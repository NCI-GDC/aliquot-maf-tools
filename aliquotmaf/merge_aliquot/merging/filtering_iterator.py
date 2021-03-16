"""
Class for filtering peekable iterator for reading in per-caller MAFs
and filtering records before even considering for overlap comparisons. This class
should be passed to the ``maflib.overlap_iter.LocatableOverlapIterator`` class.
"""


class FilteringPeekableIterator:
    """An iterator that has a `peek()` method.
    All filtering occurs on the __update_peek method.
    """

    def __init__(self, _iter):
        self._iter = _iter
        self.__update_peek()

    def __iter__(self):
        return self

    def next(self):
        """Gets the next record"""
        return self.__next__()

    def __next__(self):
        if self._peek is None:
            raise StopIteration
        to_return = self._peek
        self.__update_peek()
        return to_return

    def __update_peek(self):
        """
        Filters any variant that:
          - Isn't Somatic
          - Is larger than 50bp
          - Failed caller filters
            - Tier1,2,3,4; panel_of_normals allowed
        """
        vflt_ignore = set(["Tier1", "Tier2", "Tier3", "Tier4", "panel_of_normals"])
        can_skip = True
        while can_skip:
            self._peek = next(self._iter, None)
            if self._peek is not None:
                status = self._peek["Mutation_Status"].value.value

                size = (
                    len(self._peek["Allele"].value)
                    if self._peek["Variant_Type"].value.value == "INS"
                    else self._peek["End_Position"].value
                    - self._peek["Start_Position"].value
                )

                if status == "Somatic" and size <= 50:
                    vflt = (
                        set(
                            [i for i in self._peek["FILTER"].value if i and i != "PASS"]
                        )
                        - vflt_ignore
                    )

                    if not vflt:
                        can_skip = False
            else:
                can_skip = False

    def peek(self):
        """Returns the next element without consuming it, or None
        if there are no more elements."""
        return self._peek
