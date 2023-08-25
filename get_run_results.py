#! /usr/bin/env python3

import argparse
import boto3
import json
import logging
from datetime import datetime


parser = argparse.ArgumentParser(description='Retrieve Omics test run results')
parser.add_argument('-p', '--profile', dest='profile', metavar='PROFILE', help='AWS profile to use')
parser.add_argument('-n', '--name', dest='name', metavar='NAME', help='Name of test (run prefix)')
parser.add_argument('-v', '--verbose', dest='verbose', action='count', default=0, help='Verbosity (up to 3x)')
args = parser.parse_args()

FORMAT = '%(asctime)s %(message)s'
levels = [logging.WARNING, logging.INFO, logging.DEBUG]
level = levels[min(args.verbose, len(levels) - 1)]  # cap to last level index
logging.basicConfig(format=FORMAT, level=level)
logger = logging.getLogger('myLogger')

if (logging.root.level < logging.WARNING):
    print(f'Logging at level {logging.getLevelName(logging.root.level)}')

session = boto3.session.Session(profile_name=args.profile)
omics = session.client('omics')
s3 = session.client('s3')

# Get runs.  Note: not paginated.
ret = omics.list_runs()
all_runs = ret['items']
print(f"{len(all_runs)} runs")
# print(all_runs)
runs = [run for run in all_runs if run['name'].startswith(args.name)]
print(f"{len(runs)} runs")

# Print run summary
run0 = omics.get_run(id=runs[0]['id'])
params = run0['parameters']
print(params)

for run in runs:
    all_tasks = {}
    print(f'run {run["id"]}: {run["name"]}')
    ret = omics.get_run(id=run['id'])
    run_tasks = omics.list_run_tasks(id=run['id'])
    completed_tasks = [task for task in run_tasks['items'] if 'COMPLETED' in task['status']]
    for task in completed_tasks:
        keys = ['cpus', 'memory', 'startTime', 'stopTime']
        values = list(map(task.get, keys))
        # all_tasks[task['name']] = values
        all_tasks[task['name']] = str(values[3] - values[2])

    print(all_tasks)    