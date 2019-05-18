from OrthoEvol.Manager.database_management import DatabaseManagement


class DatabaseDispatcher(DatabaseManagement):

    def __init__(self, config_file, proj_mana, **kwargs):
        """
        The database dispatcher takes a configuration file that uses a nested dictionary of strategies that have values
        corresponding to the parameter used in the functions.  The database dispatcher inherits from the
        DatabaseManagement class so it only dispatches functions/methods from this class.

        See https://www.oreilly.com/library/view/python-cookbook/0596001673/ch01s07.html for more details.

        :param config_file:  A configuration (YAML) file that has a "Database_config" key/section.
        :type config_file:  str.
        :param proj_mana: A configuration variable for connecting projects.
        :type proj_mana: ProjectManagement.
        :param kwargs:  Key-word arguments.
        :type kwargs:  dict.
        """
        super().__init__(config_file=config_file, proj_mana=proj_mana, **kwargs)
        self.dispatcher, self.configuration = self.get_strategy_dispatcher(db_config_strategy=self.db_config_strategy)
        self.strategies = list(self.dispatcher.keys())
        self.actions = ["archive", "upload", "configure", "delete"]

    def dispatch(self, strategies=None, dispatcher=None, configuration=None):
        """
        The dispatch method takes a list of strategies and dispatches the functions using a configuration.  A nested
        dispatcher triggers the recursiveness of the dispatch function.  Dispatching means that the functions are
        serialized in an ordered dictionary, and then executed in that order.

        :param strategies:  Keys used to access the various dispatch functions and their corresponding configurations.
        :type strategies:  list.
        :param dispatcher:  A dictionary of functions that will be "dispatched".  This can be nested if using a general
        strategy.
        :type dispatcher:  dict.
        :param configuration:  A dictionary of kwargs (keyword arguments) used to configure the dispatched functions.
        :type configuration:  dict.
        """
        if not strategies:
            strategies = self.strategies
        if not dispatcher:
            dispatcher = self.dispatcher
        if not configuration:
            configuration = self.configuration

        for strategy in strategies:
            disp = dispatcher[strategy]
            conf = configuration[strategy]
            if isinstance(disp, list):
                for funk, kw in zip(disp, conf):
                    funk(**kw)
            elif isinstance(disp, dict):
                strats = list(disp.keys())
                self.dispatch(strats, disp, conf)
