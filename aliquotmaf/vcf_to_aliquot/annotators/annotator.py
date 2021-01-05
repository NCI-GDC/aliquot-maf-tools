"""
Base class for all annotators.
"""
from abc import ABCMeta, abstractmethod


class Annotator(metaclass=ABCMeta):
    def __init__(self, name=None, source=None, scheme=None):
        self.name = None
        self.source = source
        self.scheme = scheme

    @classmethod
    @abstractmethod
    def setup(cls):
        """
        Sets up and initializes the annotator instance.
        """

    @abstractmethod
    def annotate(self, **kwargs):
        """
        Performs the annotation.
        """

    @abstractmethod
    def shutdown(self):
        """
        Cleans up any connections at end of usage
        """
