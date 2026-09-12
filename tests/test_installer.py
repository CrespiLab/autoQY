import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class InstallerTests(unittest.TestCase):
    def test_batch_installer_prompts_for_and_passes_a_branch(self):
        batch = (ROOT / "Install-AutoQY.bat").read_text(encoding="utf-8")
        self.assertIn("Git branch to install [main]", batch)
        self.assertIn('-Branch "%AUTOQY_BRANCH%"', batch)

    def test_powershell_installer_validates_and_checks_remote_branch(self):
        installer = (ROOT / "Install-AutoQY.ps1").read_text(encoding="utf-8")
        self.assertIn("Select-RepositoryBranch", installer)
        self.assertIn("Assert-RemoteBranch", installer)
        self.assertIn('"ls-remote", "--exit-code", "--heads"', installer)


if __name__ == "__main__":
    unittest.main()
