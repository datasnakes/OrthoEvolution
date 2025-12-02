"""Pipeline module for parallel orthology inference workflows.

This module provides pipeline tasks that can be executed in parallel
on cluster computing systems using Luigi and SunGrid Engine.
"""

from .blastpipeline import BlastPipelineTask
from .testpipelinetask import TestPipelineTask

__all__ = ('BlastPipelineTask', 'TestPipelineTask')
