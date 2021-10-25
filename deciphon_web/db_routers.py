class DeciphonRouter:
    """
    Route deciphon unmanaged models to the deciphon database.
    (Other Django models and app-created tables use default managed database.)
    """

    deciphon_app_labels = {
        "deciphon",
    }

    def db_for_read(self, model, **hints):
        if model._meta.app_label in self.deciphon_app_labels:
            return "deciphon"
        return None

    def db_for_write(self, model, **hints):
        if model._meta.app_label in self.deciphon_app_labels:
            return "deciphon"
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Explicitly allow relations between deciphon models. Otherwise let Django figure it out.
        Cross-database relations are generally not supported.
        """
        if (
            obj1._meta.app_label in self.deciphon_app_labels
            or obj2._meta.app_label in self.deciphon_app_labels
        ):
            return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        No migrations for deciphon models – they are externally managed.
        """
        if app_label in self.deciphon_app_labels:
            return False
        return None
