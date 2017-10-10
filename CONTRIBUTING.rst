====================================
Contributing to this project
====================================

Thanks for your desire to contribute to this project.

Preparing your Fork
^^^^^^^^^^^^^^^^^^^
1. Hit 'fork' on Github, creating e.g. ``yourname/Datasnakes-Scripts``.
2. Clone your project: ``git clone git@github.com:yourname/Datasnakes-Scripts``.
3. Create a branch: ``cd Datasnakes-Scripts; git checkout -b foo-the-bars 1.3``.

Development Mode Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Change to the project repository on your machine: ``cd Datasnakes-Scripts``
2. Install in development mode: ``pip install -e .``


Making your Changes
^^^^^^^^^^^^^^^^^^^
1. Add your contributions
2. Run tests and make sure they pass. (We use unittest.)
3. Commit your changes: ``git commit -m "Foo the bars"``



Creating Pull Requests
^^^^^^^^^^^^^^^^^^^^^^

1. Push your commit to get it back up to your fork: ``git push origin HEAD``
2. Visit Github, click handy "Pull request" button that it will make upon
   noticing your new branch.
3. In the description field, write down issue number (if submitting code fixing
   an existing issue) or describe the issue + your fix (if submitting a wholly
   new bugfix).
4. Hit 'submit'! And please be patient - we will get to you when
   we can.

Version control branching
-------------------------

* Always **make a new branch** for your work, no matter how small. This makes
  it easy for others to take just that one set of changes from your repository,
  in case you have multiple unrelated changes floating around.

    * A corollary: **don't submit unrelated changes in the same branch/pull
      request**! The maintainer shouldn't have to reject your awesome bugfix
      because the feature you put in with it needs more review.

* **Base your new branch off of the appropriate branch** on the main
  repository:

    * **Bug fixes** should be based on the branch named after the **oldest
      supported release line** the bug affects.

        * E.g. if a feature was introduced in 1.1, the latest release line is
          1.3, and a bug is found in that feature - make your branch based on
          1.1.  The maintainer will then forward-port it to 1.3 and master.
        * Bug fixes requiring large changes to the code or which have a chance
          of being otherwise disruptive, may need to base off of **master**
          instead. This is a judgement call -- ask the devs!

    * **New features** should branch off of **the 'master' branch**.

        * Note that depending on how long it takes for the dev team to merge
          your patch, the copy of ``master`` you worked off of may get out of
          date! If you find yourself 'bumping' a pull request that's been
          sidelined for a while, **make sure you rebase or merge to latest
          master** to ensure a speedier resolution.

Code formatting
---------------

* **Follow the style you see used in the primary repository**! Consistency with
  the rest of the project always trumps other considerations.
* Python projects usually follow the `PEP-8
  <http://www.python.org/dev/peps/pep-0008/>`_ guidelines (though many have
  minor deviations depending on the lead maintainers' preferences.)

Documentation isn't optional
----------------------------

It's not! Patches without documentation will be returned to sender.  By
"documentation" we mean:

* **Docstrings** (for Python; or API-doc-friendly comments for other languages)
  must be created or updated for public API functions/methods/etc. (This step
  is optional for some bugfixes.)
* New features should ideally include updates to **prose documentation**,
  including useful example code snippets.

Tests aren't optional
---------------------

Any bugfix that doesn't include a test proving the existence of the bug being
fixed, may be suspect.  Ditto for new features that can't prove they actually
work.

We've found that test-first development really helps make features better
architected and identifies potential edge cases earlier instead of later.
Writing tests before the implementation is strongly encouraged.



Thanks to https://contribution-guide-org.readthedocs.io/ for the general guide.