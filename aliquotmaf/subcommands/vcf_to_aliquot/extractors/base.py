"""Abstract base class for extracting data"""
from abc import abstractmethod
from typing import Any, Protocol


class Extractor(Protocol):
    @classmethod
    @abstractmethod
    def extract(cls, *args: Any, **kwargs: Any) -> Any:
        """
        All extractors much implement the extraction function.
        """
