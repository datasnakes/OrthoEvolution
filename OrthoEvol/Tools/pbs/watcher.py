import os
import asyncio
from pathlib import Path
from time import sleep, time
from datetime import datetime
import importlib.util
import argparse
import yaml
import qwatch

parser = argparse.ArgumentParser(description='A job watcher.')

parser.add_argument('-yamlfile', action="store", dest="yamlfile")
_ = parser.parse_args()
yfile = _.yamlfile
with open(yfile, 'r') as yf:
    _kwargs = yaml.load(yf)


async def _async_watch(job_id, directory, python_datetime=None, first_time=True, **kwargs):
    """Wait until a list of jobs finishes and get updates."""

    watch_one = qwatch.Qwatch(jobs=[job_id], directory=directory,  watch=True, users=None, **kwargs)
    job_dict = watch_one.update_qstat_data(save=True, process=True, data=True)

    if job_dict:
        data = watch_one.get_data(data_frame=True, python_datetime=python_datetime)
        watch_one.update_csv(file=watch_one.data_filename, data=data)
        print(f"Updated qstat data for {job_id}")

        if job_dict[job_id]['job_state'] == 'Q':
            await asyncio.sleep(watch_one.sleeper)
            await _async_watch(job_id, directory, datetime.now(), first_time=True **kwargs)
        elif job_dict[job_id]['job_state'] == 'R':
            # Create the static data file on the first instance of a running job
            if first_time:
                info = watch_one.get_info(data_frame=True, python_datetime=python_datetime)
                watch_one.update_csv(file=watch_one.info_filename, data=info)
            await asyncio.sleep(watch_one.sleeper)
            await _async_watch(job_id, directory, datetime.now(), first_time=False, **kwargs)
    watch_one.plot_memory()

    return f'Finished {job_id}'


def get_watcher_tasks(jobs, home_dir, **kwargs):
    python_datetime = datetime.now()
    tasks = [asyncio.ensure_future(_async_watch(job_id=job, directory=home_dir/Path(job),
                                                python_datetime=python_datetime, **kwargs)) for job in jobs]
    return tasks


ioloop = asyncio.get_event_loop()
tasks = get_watcher_tasks(jobs=_kwargs["jobs"],
                          home_dir=_kwargs["directory"], **_kwargs["kwargs"])
ioloop.run_until_complete(asyncio.wait(tasks))
ioloop.close()
