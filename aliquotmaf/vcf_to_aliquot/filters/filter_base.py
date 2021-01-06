#!/usr/bin/env python3
"""
Base class for all filters.
"""
from abc import ABCMeta, abstractmethod


class Filter(metaclass=ABCMeta):
    def __init__(self, name=None, source=None):
        self.name = None
        self.source = source
        self.tags = []

    @classmethod
    @abstractmethod
    def setup(cls):
        """
        Sets up and initializes the filter instance.
        """

    @abstractmethod
    def filter(self, **kwargs):
        """
        Performs the filter.
        """

    @abstractmethod
    def shutdown(self):
        """
        Cleans up any connections at end of usage
        """


# __END__
