import unittest
import base64
import json
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).parents[1]


class InstallerTests(unittest.TestCase):
    def test_batch_installer_prompts_for_and_passes_a_branch(self):
        batch = (ROOT / "Install-AutoQY.bat").read_text(encoding="utf-8")
        self.assertIn("Git branch to install [main]", batch)
        self.assertIn('-Branch "%AUTOQY_BRANCH%"', batch)
        self.assertIn(
            "raw.githubusercontent.com/CrespiLab/autoQY/%AUTOQY_BRANCH%/Install-AutoQY.ps1",
            batch,
        )
        self.assertNotIn(
            "raw.githubusercontent.com/CrespiLab/autoQY/main/Install-AutoQY.ps1",
            batch,
        )
        self.assertIn('-InstallerSourceUrl "%AUTOQY_PS1_URL%"', batch)
        self.assertIn("copy and paste its full path", batch)
        self.assertIn("Press Enter to close", batch)

    def test_powershell_installer_validates_and_checks_remote_branch(self):
        installer = (ROOT / "Install-AutoQY.ps1").read_text(encoding="utf-8")
        self.assertIn("Select-RepositoryBranch", installer)
        self.assertIn("Assert-RemoteBranch", installer)
        self.assertIn('"ls-remote", "--exit-code", "--heads"', installer)
        self.assertIn("Get-InstallerSourceUrl", installer)
        self.assertIn("PowerShell installer URL", installer)
        self.assertIn("Enter or paste the full path", installer)
        self.assertIn('Read-Host "Press Enter to close"', installer)


@unittest.skipUnless(os.name == "nt" and shutil.which("powershell.exe"),
                     "Windows PowerShell is required for launcher integration tests")
class InstallerLauncherTests(unittest.TestCase):
    def run_launcher(self, folder, installer_exit=0, download_failure=False):
        """Run the real BAT flow, replacing only its network download with a stub."""
        batch = (ROOT / "Install-AutoQY.bat").read_text(encoding="utf-8")
        installer = """param([string]$Branch, [string]$InstallerSourceUrl, [switch]$NoClosePrompt)
$details = @{branch=$Branch; source=$InstallerSourceUrl; cwd=(Get-Location).Path}
[IO.File]::WriteAllText((Join-Path (Get-Location) 'invocation.json'), ($details | ConvertTo-Json))
exit INSTALLER_EXIT
# PADDING
""".replace("INSTALLER_EXIT", str(installer_exit)).replace("PADDING", "x" * 1100)
        payload = base64.b64encode(installer.encode("utf-8")).decode("ascii")
        mock_download = (
            "function Invoke-WebRequest { param($Uri, $OutFile, [switch]$UseBasicParsing); "
            f"[IO.File]::WriteAllBytes($OutFile, [Convert]::FromBase64String('{payload}')); "
            + ("throw 'Simulated download failure'; " if download_failure else "") + "}; "
        )
        lines = batch.splitlines()
        download_line = next(index for index, line in enumerate(lines)
                             if "Invoke-WebRequest -Uri" in line)
        lines[download_line] = lines[download_line].replace(
            '-Command "', '-Command "' + mock_download, 1)
        launcher = folder / "Install-AutoQY.bat"
        launcher.write_text("\n".join(lines) + "\n", encoding="utf-8")
        existing = folder / "Install-AutoQY.ps1"
        existing.write_text("Original repository installer; preserve me.", encoding="utf-8")
        temporary = folder / "temp files"
        temporary.mkdir()
        environment = dict(os.environ, TEMP=str(temporary), TMP=str(temporary))
        result = subprocess.run(
            ["cmd.exe", "/d", "/c", str(launcher)], input="feature/test\n\n",
            cwd=folder, env=environment, capture_output=True, text=True, timeout=45,
        )
        self.assertEqual(existing.read_text(), "Original repository installer; preserve me.")
        self.assertEqual(list(temporary.iterdir()), [], result.stdout + result.stderr)
        return result

    def test_success_preserves_existing_installer_and_cleans_only_temporary_download(self):
        with TemporaryDirectory(prefix="AutoQY installer test ") as directory:
            folder = Path(directory)
            result = self.run_launcher(folder)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            details = json.loads((folder / "invocation.json").read_text())
            self.assertEqual(details["branch"], "feature/test")
            self.assertIn("/feature/test/Install-AutoQY.ps1", details["source"])
            self.assertEqual(Path(details["cwd"]), folder)

    def test_installer_failure_preserves_exit_code_and_removes_download(self):
        with TemporaryDirectory(prefix="AutoQY installer test ") as directory:
            result = self.run_launcher(Path(directory), installer_exit=7)
            self.assertEqual(result.returncode, 7, result.stdout + result.stderr)

    def test_failed_download_cleans_partial_file_and_never_runs_installer(self):
        with TemporaryDirectory(prefix="AutoQY installer test ") as directory:
            folder = Path(directory)
            result = self.run_launcher(folder, download_failure=True)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertFalse((folder / "invocation.json").exists())


if __name__ == "__main__":
    unittest.main()
