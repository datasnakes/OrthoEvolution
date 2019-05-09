from OrthoEvol.Manager.database_management import DatabaseManagement


class DatabaseDispatcher(DatabaseManagement):

    def __init__(self, config_file, proj_mana, upload_refseq_release=False, **kwargs):
        super().__init__(config_file=config_file, proj_mana=proj_mana)
        self.dispatcher, self.configuration = self.get_strategy_dispatcher(db_config_strategy=self.db_config_strategy)
        self.strategies = list(self.dispatcher.keys())
        self.actions = ["archive", "upload", "configure", "delete"]

        if upload_refseq_release == True:
            self.refseq_release_dispatcher(**kwargs)

    def dispatch(self, strategies, dispatcher, configuration):
        for strategy in strategies:
            disp = dispatcher[strategy]
            conf = configuration[strategy]
            if isinstance(disp, list):
                for funk, kw in zip(disp, conf):
                    funk(**kw)
            elif isinstance(disp, dict):
                strats = list(disp.keys())
                self.dispatch(strats, disp, conf)

    def uploading_dispatch(self, strategies):
        dispatcher = {}
        configuration = {}
        for strategy in strategies:
            dispatcher[strategy] = []
            configuration[strategy] = []
            if strategy == "NCBI_refseq_release":
                self.configuration[strategy]['configure_flag'] = False
                self.configuration[strategy]['upload_flag'] = True
                disp, conf = self.NCBI_refseq_release(**self.configuration[strategy])
                dispatcher[strategy].append(disp)
                configuration[strategy].append(conf)

        self.dispatch(strategies, dispatcher, configuration)

    def refseq_release_dispatcher(self, **kwargs):
        self.upload_refseq_release_files(**kwargs)


