import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
