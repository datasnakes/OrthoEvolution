from OrthoEvol.Manager.database_management import DatabaseManagement


class DatabaseDispatcher(DatabaseManagement):

    def __init__(self, config_file, proj_mana, upload_refseq_release=False, **kwargs):
        super().__init__(config_file=config_file, proj_mana=proj_mana)
        self.dispatcher, self.configuration = self.get_strategy_dispatcher(db_config_strategy=self.db_config_strategy)
        self.strategies = list(self.dispatcher.keys())
        self.actions = ["archive", "upload", "configure", "delete"]
        # if len(self.dispatcher.keys()) == 1:
        #     self.strategy = list(self.dispatcher.keys())[0]
        #
        #     # Create class methods from dispatcher
        #     self.archive_disp = self.dispatcher[self.strategy]["archive"]
        #     self.configure_disp = self.dispatcher[self.strategy]["configure"]
        #     self.upload_disp = self.dispatcher[self.strategy]["upload"]
        #
        #     # Create class attributes from configuration
        #     self.archive_config = self.configuration[self.strategy]["archive"]
        #     self.configure_config = self.configuration[self.strategy]["configure"]
        #     self.upload_config = self.configuration[self.strategy]["configure"]
        # else:
        #     for disp, c in zip(self.dispatcher, self.configuration):

        if upload_refseq_release == True:
            self.refseq_release_dispatcher(**kwargs)

    def dispatch(self, strategy, dispatcher, configuration):
        disp = dispatcher[strategy]
        conf = configuration[strategy]
        file_list = None
        if isinstance(disp, list):
            for funk, kw in zip(disp, conf):
                files = funk(file_list, **kw)
                if isinstance(files, list):
                    file_list = files
        elif isinstance(disp, dict):
            for action in disp.keys():
                self.dispatch(action, disp, conf)

    def refseq_release_dispatcher(self, **kwargs):
        self.upload_refseq_release_files(kwargs)


