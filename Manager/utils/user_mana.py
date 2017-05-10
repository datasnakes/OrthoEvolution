import os
from pathlib import Path

from Manager.utils.repo_mana import RepoMana as RM


class UserMana(RM):

    def __init__(self, user, home=os.getcwd(), project=None, new_project=False):
        super().__init__(home=home, user=user)

