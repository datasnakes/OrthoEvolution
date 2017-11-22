from OrthoEvol.Manager.database_management import DatabaseManagement


class DatabaseDispatcher(DatabaseManagement):

    def __init__(self, config_file, **kwargs):
        super().__init__(config_file=config_file)
        self.dispatcher, self.configuration = self.get_strategy_dispatcher(db_config_strategy=kwargs)
        if len(self.dispatcher.keys()) == 1:
            self.strategy = list(self.dispatcher.keys())[0]
        else:
            raise ValueError

        self.archive_disp = self.dispatcher[self.strategy]["archive"]
        self.configure_disp = self.dispatcher[self.strategy]["configure"]
        self.upload_disp = self.dispatcher[self.strategy]["upload"]

        self.archive_config = self.configuration[self.strategy]["archive"]
        self.configure_config = self.configuration[self.strategy]["configure"]
        self.upload_config = self.configuration[self.strategy]["configure"]

    def dispatch_the_uploader(self):
        self.upload_disp(**self.upload_config)


