#!/usr/bin/python3
from typing import List

class Partitioner:
    def __init__(
        self,
        script: str,
        format: str | List[str],
        *,
        parallel: bool,
        dynamic_header: bool = False,
    ) -> None:
        self.script = script
        if isinstance(format, str):
            self.format = [format]
        else:
            self.format = format
        self.parallel = parallel
        self.dynamic_header = dynamic_header


partitioner_mapping = {
    "Mt-KaHyPar":       Partitioner("mt_kahypar", ["graph", "hmetis"], parallel=True, dynamic_header=True)
}