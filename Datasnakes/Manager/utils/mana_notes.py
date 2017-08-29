<<<<<<< HEAD
import os
from Datasnakes.Manager.utils import ProjMana, WebMana, UserMana, RepoMana, Mana, DataMana
repo = "Test"
user = "rgilmore"
project = "Orthlogs_Test"
research = "comparative genetics"
research_type = "public"
website = 'fun'




ProjMana(repo=repo, user=user, project=project, research=None, research_type=None, app=None, home=os.getcwd(), new_project=False, new_research=False, new_app=False, **kwargs)

WebMana(repo=repo, website=website, host='0.0.0.0', port='5252', home=os.getcwd(), new_website=False, create_admin=False, **kwargs)

UserMana(repo=repo, user=user, project=None, database=None, home=os.getcwd(), new_user=False, new_project=False, new_db=False, **kwargs)

RepoMana(repo=repo, user=None, home=os.getcwd(), new_user=False, new_repo=False)

Mana(repo=None, home=os.getcwd(), new_repo=False)


DataMana(home, data, new_data, repo, new_repo, user, new_user, project, new_project, app, new_app, database, new_db,
         website, new_website, host, port, create_admin, research, new_research, research_type)

=======
>>>>>>> b4e6bb4ddfa7bc087f1fae5a8844594a6a6198c4
