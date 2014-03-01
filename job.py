

class Job(object):

    def __init__(self, **kwargs):

        for key in kwargs["code"].items():
            setattr(self, key, val)

