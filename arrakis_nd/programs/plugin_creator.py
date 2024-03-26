"""
Script for generating a plugin template with an accompanying config file
"""
import os
import shutil
import argparse
import socket

from arrakis_nd.utils.utils import get_datetime


def run():
    parser = argparse.ArgumentParser(
        prog='Arrakis Plugin Creator',
        description='This program constructs a Arrakis Plugin ' +
                    'with an empty template.',
        epilog='...'
    )
    parser.add_argument(
        '-plugin_name', dest='plugin_name', default='my',
        help='name to use for the custom plugin parts.'
    )
    parser.add_argument(
        '-plugin_location', dest='plugin_location', default='/local_scratch',
        help='location for the local scratch directory.'
    )

    args = parser.parse_args()
    if args.plugin_location is not None:
        if not os.path.isdir(args.plugin_location):
            args.plugin_location = './'
    else:
        args.plugin_location = './'
    snake_case_name = args.plugin_name
    plugin_name_words = snake_case_name.split('_')
    camel_case_name = ''.join(word.capitalize() for word in plugin_name_words)
    plugin_location = args.plugin_location + '/' + camel_case_name + 'CustomPlugin/'

    # now copy template files and replace 'Empty' with args.plugin_name
    now = get_datetime()
    if not os.path.isdir(plugin_location):
        os.makedirs(plugin_location)
    else:
        if not os.path.isdir(f"{plugin_location}.backup/"):
            os.makedirs(f"{plugin_location}.backup/")
        os.makedirs(f"{plugin_location}.backup/{now}")
        selected_files = [
            file for file in os.listdir(plugin_location)
            if (file.endswith('.py') or file.endswith('.yaml'))
        ]
        for file in selected_files:
            source_path = os.path.join(plugin_location, file)
            destination_path = os.path.join(f"{plugin_location}.backup/{now}", file)
            shutil.move(source_path, destination_path)

    # add *user_name* and *user_date*
    try:
        user_name = socket.gethostname()
    except Exception:
        user_name = "*user_name*"

    shutil.copy(
        os.path.dirname(__file__) + '/../plugins/empty_plugin.py',
        f'{plugin_location}/{snake_case_name}.py'
    )
    with open(f'{plugin_location}/{snake_case_name}.py', 'r') as file:
        content = file.read()
    modified_content = content.replace(
        'Empty',
        camel_case_name
    )
    modified_content = modified_content.replace(
        'empty',
        snake_case_name
    )
    modified_content = modified_content.replace(
        "*user_name*", user_name
    )
    modified_content = modified_content.replace(
        "*user_date*", now
    )
    with open(f'{plugin_location}{snake_case_name}.py', 'w') as file:
        file.write(modified_content)


if __name__ == "__main__":
    run()
