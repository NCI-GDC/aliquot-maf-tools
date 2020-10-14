"""Abstract base class for extracting data"""
from abc import ABCMeta, abstractmethod


class Extractor(metaclass=ABCMeta):
    @classmethod
    @abstractmethod
    def extract(cls, **kwargs):
        """
        All extractors much implement the extraction function.
        """
