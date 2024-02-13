class Checker:
    """Utility to hold checks


    The method :meth:`~sequana.utils.checks.Checke./tryme` calls the method of function
    provide. This method is expected to return a dictionary
    with 2 keys called status and msg. Status should be
    in 'Error', 'Warning', 'Success'.

    The attributes hold all calls to :meth:`tryme`


    """

    def __init__(self):

        self.results = []

    def tryme(self, meth):
        try:
            status = meth()
            if "msg" in status and "status" in status:
                self.results.append(status)
            else:
                self.results.append({"status": "Success", "msg": status})
        except Exception as err:  # pragma: no cover
            self.results.append({"msg": err, "status": "Error"})
