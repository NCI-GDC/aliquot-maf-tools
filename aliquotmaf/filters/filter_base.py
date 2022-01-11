"""
Base class for all filters.
"""
import logging
from abc import ABC, abstractmethod
from typing import Any, List, Optional

from aliquotmaf.logger import Logger


class Filter(ABC):
    name: str
    source: Optional[str]
    logger: logging.Logger
    tags: List[str] = []

    def __init__(self, name: str, source: Optional[str] = None):
        self.name = name
        self.source = source
        self.logger = Logger.get_logger(self.__class__.__name__)
        self.tags: List[str] = []

    @classmethod
    @abstractmethod
    def setup(cls, *args: Any, **kwargs: Any) -> 'Filter':
        """
        Sets up and initializes the filter instance.
        """

    @abstractmethod
    def filter(self, *args: Any, **kwargs: Any) -> bool:
        """
        Performs the filter.
        """

    @abstractmethod
    def shutdown(self) -> None:
        """
        Cleans up any connections at end of usage
        """
