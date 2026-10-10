"""Run Windows uninstall checks without deleting files or killing real processes."""
import json
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).parents[1]


def ps_string(value):
    return "'" + str(value).replace("'", "''") + "'"


class UninstallerSourceTests(unittest.TestCase):
    def test_new_shortcuts_record_the_exact_environment(self):
        installer = (ROOT / "Install-AutoQY.ps1").read_text(encoding="utf-8")
        self.assertIn('-EnvironmentPath `"$EnvironmentPath`"', installer)


@unittest.skipUnless(os.name == "nt" and shutil.which("powershell.exe"),
                     "Windows PowerShell is required for uninstall integration checks")
class UninstallerWindowsTests(unittest.TestCase):
    def run_powershell(self, program, directory):
        wrapper = directory / "check.ps1"
        wrapper.write_text(program, encoding="utf-8-sig")
        return subprocess.run(
            ["powershell.exe", "-NoLogo", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-File", str(wrapper)], capture_output=True, text=True, timeout=30,
        )

    def test_process_checks_are_scoped_and_do_not_kill_anything(self):
        program = r"""
$ErrorActionPreference = 'Stop'
$parseTokens = $null
$parseErrors = $null
$ast = [Management.Automation.Language.Parser]::ParseFile(SCRIPT, [ref]$parseTokens, [ref]$parseErrors)
if ($parseErrors.Count) { throw ($parseErrors -join '\n') }
foreach ($definition in $ast.FindAll({ param($node)
    $node -is [Management.Automation.Language.FunctionDefinitionAst]
}, $false)) { Invoke-Expression $definition.Extent.Text }
$script:Processes = @(
    [pscustomobject]@{ProcessId=$PID; Name='powershell.exe'; ExecutablePath='C:\Windows\powershell.exe'; CommandLine='uninstall C:\AutoQY-Core'},
    [pscustomobject]@{ProcessId=11; Name='python.exe'; ExecutablePath='C:\other\python.exe'; CommandLine='python C:\AutoQY-Core\server.py'},
    [pscustomobject]@{ProcessId=12; Name='python.exe'; ExecutablePath='C:\other\python.exe'; CommandLine='python C:\AutoQY-Core-copy\server.py'},
    [pscustomobject]@{ProcessId=13; Name='pythonw.exe'; ExecutablePath='C:\envs\autoqy\pythonw.exe'; CommandLine='pythonw -m autoqy_core.tools.smoother_gui'},
    [pscustomobject]@{ProcessId=14; Name='python.exe'; ExecutablePath='C:\envs\autoqy\python.exe'; CommandLine='python unrelated.py'},
    [pscustomobject]@{ProcessId=15; Name='autoqy-analysis-gui.exe'; ExecutablePath='C:\envs\different\Scripts\autoqy-analysis-gui.exe'; CommandLine='autoqy-analysis-gui.exe'},
    [pscustomobject]@{ProcessId=16; Name='autoqy-power-gui.exe'; ExecutablePath=$null; CommandLine=$null},
    [pscustomobject]@{ProcessId=17; Name='autoqy-smoother-gui.exe'; ExecutablePath='C:\envs\autoqy-copy\Scripts\autoqy-smoother-gui.exe'; CommandLine='autoqy-smoother-gui.exe'}
)
function Get-CimInstance { param($ClassName, $ErrorAction); $script:Processes }
function Stop-Process { throw 'Must not kill user processes' }
$keep = @(Get-UninstallBlockingProcesses 'C:\AutoQY-Core' 'C:\envs\autoqy' $false | ForEach-Object ProcessId)
$remove = @(Get-UninstallBlockingProcesses 'C:\AutoQY-Core' 'C:\envs\autoqy' $true | ForEach-Object ProcessId)
@{keep=$keep; remove=$remove} | ConvertTo-Json -Compress
""".replace("SCRIPT", ps_string(ROOT / "Uninstall-AutoQY.ps1"))
        with TemporaryDirectory(prefix="AutoQY uninstall check ") as temporary:
            result = self.run_powershell(program, Path(temporary))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        actual = json.loads(result.stdout)
        self.assertEqual(actual["keep"], [11, 13, 16])
        self.assertEqual(actual["remove"], [11, 13, 14, 16])

    def exercise_uninstaller(self, *, active=False, cancel=False, conda_exit=0,
                             scan_failure=False, old_shortcut=False, keep_environment=False,
                             environment_absent=False, reopened=False):
        with TemporaryDirectory(prefix="AutoQY uninstall flow ") as temporary:
            folder = Path(temporary)
            installation = folder / "AutoQY-Core"
            installation.mkdir()
            (installation / "pyproject.toml").write_text('[project]\nname="autoqy-core"\n')
            (installation / "autoqy_core").mkdir()
            environment = folder / "envs" / "autoqy-test"
            if not environment_absent:
                environment.mkdir(parents=True)
            event_file = folder / "events.json"
            program = r"""
$ErrorActionPreference = 'Stop'
$global:UninstallTestEvents = [Collections.Generic.List[string]]::new()
$global:UninstallTestActive = ACTIVE
function Read-Host {
    param($Prompt)
    if ($Prompt -like 'Are you sure*') { return 'y' }
    if ($Prompt -like 'Also remove*') { return KEEP_ANSWER }
    if ($Prompt -like 'Press Enter to check*') {
        if (CANCEL) { $global:UninstallTestEvents.Add('cancel'); return 'cancel' }
        $global:UninstallTestEvents.Add('closed-gui')
        $global:UninstallTestActive = $false
    }
    return ''
}
function Get-CimInstance {
    param($ClassName, $ErrorAction)
    $global:UninstallTestEvents.Add('process-check')
    if (SCAN_FAILURE) { throw 'Simulated Windows process-check failure' }
    if ($global:UninstallTestActive) {
        [pscustomobject]@{ProcessId=101; Name='pythonw.exe'; ExecutablePath=PYTHON; CommandLine='pythonw -m autoqy_core.tools.smoother_gui'}
    }
}
function Mock-Conda {
    param([string[]]$Arguments)
    if ($Arguments[1] -eq 'list') {
        $global:UninstallTestEvents.Add('environment-list')
        $global:LASTEXITCODE = 0
        @{envs=@(ENVIRONMENT)} | ConvertTo-Json -Compress
        return
    }
    if ($global:UninstallTestActive) { throw 'Conda was called before the GUI closed' }
    $global:UninstallTestEvents.Add('environment-remove')
    if ($Arguments[2] -ne '--prefix' -or $Arguments[3] -ne ENVIRONMENT) {
        throw 'Conda removal must target the exact recorded environment'
    }
    $global:LASTEXITCODE = CONDA_EXIT
    if (REOPENED) {
        $global:UninstallTestEvents.Add('reopened-gui')
        $global:UninstallTestActive = $true
    }
}
function Remove-Item {
    param($LiteralPath, [switch]$Recurse, [switch]$Force)
    if ($global:UninstallTestActive) { throw 'Files were removed before the GUI closed' }
    $global:UninstallTestEvents.Add('remove:' + $LiteralPath)
}
try {
    & SCRIPT -ProjectRoot INSTALLATION -EnvironmentName 'autoqy-test' `
        -CondaCommand 'Mock-Conda' -ShortcutDirectory SHORTCUTS ENVIRONMENT_ARGUMENT
    $resultCode = $LASTEXITCODE
} finally {
    [IO.File]::WriteAllText(EVENT_FILE, (ConvertTo-Json -InputObject @($global:UninstallTestEvents)))
}
exit $resultCode
"""
            replacements = {
                "ACTIVE": "$true" if active else "$false",
                "CANCEL": "$true" if cancel else "$false",
                "SCAN_FAILURE": "$true" if scan_failure else "$false",
                "KEEP_ANSWER": "'n'" if keep_environment else "'y'",
                "CONDA_EXIT": str(conda_exit),
                "REOPENED": "$true" if reopened else "$false",
                "PYTHON": ps_string(environment / "pythonw.exe"),
                "SCRIPT": ps_string(ROOT / "Uninstall-AutoQY.ps1"),
                "INSTALLATION": ps_string(installation),
                "SHORTCUTS": ps_string(folder / "shortcuts"),
                "ENVIRONMENT_ARGUMENT": "" if old_shortcut else "-EnvironmentPath " + ps_string(environment),
                "ENVIRONMENT": ps_string(environment),
                "EVENT_FILE": ps_string(event_file),
            }
            for token in sorted(replacements, key=len, reverse=True):
                program = program.replace(token, replacements[token])
            result = self.run_powershell(program, folder)
            events = json.loads(event_file.read_text())
            # Remove-Item is mocked: only the disposable fixture is touched.
            self.assertTrue((installation / "pyproject.toml").exists())
            return result, events

    def test_open_gui_waits_for_closure_before_any_removal(self):
        result, events = self.exercise_uninstaller(active=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Save your work and close all AutoQY windows", result.stdout)
        self.assertLess(events.index("closed-gui"), events.index("environment-remove"))
        self.assertTrue(any(event.startswith("remove:") for event in events))
        self.assertEqual(events.count("process-check"), 3)

    def test_cancel_with_open_gui_leaves_environment_and_files_alone(self):
        result, events = self.exercise_uninstaller(active=True, cancel=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("no files or environments were changed", result.stdout)
        self.assertNotIn("environment-remove", events)
        self.assertFalse(any(event.startswith("remove:") for event in events))

    def test_conda_failure_preserves_installation_and_shortcuts(self):
        result, events = self.exercise_uninstaller(conda_exit=7)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("installation folder and shortcuts were kept", result.stdout)
        self.assertFalse(any(event.startswith("remove:") for event in events))

    def test_process_check_failure_stops_before_deletion(self):
        result, events = self.exercise_uninstaller(scan_failure=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("process-check failure", result.stdout)
        self.assertNotIn("environment-remove", events)
        self.assertFalse(any(event.startswith("remove:") for event in events))

    def test_old_shortcut_resolves_environment_and_checks_open_gui(self):
        result, events = self.exercise_uninstaller(active=True, old_shortcut=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("environment-list", events)
        self.assertLess(events.index("closed-gui"), events.index("environment-remove"))

    def test_keep_environment_still_waits_for_autoqy_to_close(self):
        result, events = self.exercise_uninstaller(active=True, keep_environment=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("closed-gui", events)
        self.assertNotIn("environment-remove", events)
        self.assertTrue(any(event.startswith("remove:") for event in events))

    def test_already_removed_environment_does_not_block_folder_cleanup(self):
        result, events = self.exercise_uninstaller(environment_absent=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("already absent", result.stdout)
        self.assertNotIn("environment-remove", events)
        self.assertTrue(any(event.startswith("remove:") for event in events))

    def test_reopened_gui_is_checked_again_before_folder_removal(self):
        result, events = self.exercise_uninstaller(reopened=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        removal = next(index for index, event in enumerate(events) if event.startswith("remove:"))
        self.assertLess(events.index("reopened-gui"), events.index("closed-gui"))
        self.assertLess(events.index("closed-gui"), removal)


if __name__ == "__main__":
    unittest.main()
