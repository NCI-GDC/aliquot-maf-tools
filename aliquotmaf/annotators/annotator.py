"""
Base class for all annotators.
"""
from abc import ABCMeta, abstractmethod
from typing import TYPE_CHECKING

from aliquotmaf.logger import Logger

if TYPE_CHECKING:
    from maflib.record import MafRecord


class Annotator(metaclass=ABCMeta):
    def __init__(self, name=None, source=None, scheme=None):
        self.name = None
        self.source = source
        self.scheme = scheme
        self.logger = Logger.get_logger(self.__class__.__name__)

    @classmethod
    @abstractmethod
    def setup(cls):
        """
        Sets up and initializes the annotator instance.
        """

    @abstractmethod
    def annotate(self, *args, **kwargs) -> 'MafRecord':
        """
        Performs the annotation.
        """

    @abstractmethod
    def shutdown(self):
        """
        Cleans up any connections at end of usage
        """
