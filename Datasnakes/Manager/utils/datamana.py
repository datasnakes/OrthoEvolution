from Manager.utils.mana import ProjMana as PM


class DataMana(PM):

    def __init__(self, repo, user, project, research, research_type, data, new_data):
        super().__init__(repo=repo, user=user, project=project, research=research, research_type=research_type)