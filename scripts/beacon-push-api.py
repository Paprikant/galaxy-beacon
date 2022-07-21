# -*- coding: utf-8 -*-
"""Pushes genomic variant information from a galaxy to a beacon instance.

Usage:
    ./beacon-push.py
"""

import argparse
import asyncio
import datetime
import json
import logging
import os
import re
from dataclasses import dataclass
from typing import Dict, List

# import utilities from beacon-python
# pip install git+https://github.com/CSCfi/beacon-python
#
# TODO can we install it like this
from beacon_api.utils.db_load import BeaconDB
from bioblend.galaxy import GalaxyInstance
# cyvcf2 is a vcf parser
from cyvcf2 import VCF, Variant


def parse_arguments():
    """
    Parses command line arguments with argparse
    """
    parser = argparse.ArgumentParser(description="Push genomic variants from galaxy to beacon.")

    # arguments controlling output 
    parser.add_argument("-v", "--verbosity", action="count", default=0,
                        help="log verbosity, can be repeated up to three times")

    # arguments controlling galaxy connection
    parser.add_argument("-u", "--galaxy-url", type=str, metavar="", default="http://localhost:8080", dest="galaxy_url",
                        help="galaxy hostname or IP")
    parser.add_argument("-k", "--galaxy-key", type=str, metavar="", default="6edbc8a89bbff89bb5232867edc1183c",
                        dest="galaxy_key", help="API key of a galaxy user WITH ADMIN PRIVILEGES")

    return parser.parse_args()


def set_up_logging(verbosity: int):
    """
    Configures the logger for this script

        Parameters:
            verbosity (int):
    """

    # configure log level to match the given verbosity
    if verbosity > 1:
        logging.basicConfig(level=logging.DEBUG)
    elif verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)


def set_up_galaxy_instance(galaxy_url: str, galaxy_key: str) -> GalaxyInstance:
    """
    Returns a galaxy instance with the given URL and api key.
    Exits immediately if either connection or authentication to galaxy fails.

        Parameters:
            galaxy_url (str): Base URL of a galaxy instance
            galaxy_key (str): API key with admin privileges

        Returns:
            gi (GalaxyInstance): Galaxy instance with confirmed admin access to the given galaxy instance
    """

    logging.info("trying to connect to galaxy at {}".format(galaxy_url))

    # configure a galaxy instance with galaxy_url and api_key
    try:
        gi = GalaxyInstance(galaxy_url, key=galaxy_key)
    except Exception as e:
        # if galaxy_url does not follow the scheme <protocol>://<host>:<port>, GalaxyInstance attempts guessing the URL 
        # this exception is thrown when neither "http://<galaxy_url>:80" nor "https://<galaxy_url>:443" are accessible
        logging.critical("failed to guess URL from \"{}\" - {}".format(galaxy_url, e))
        exit(2)

    # test network connection and successful authentification
    try:
        response = gi.make_get_request(galaxy_url + "/api/whoami")
        content = json.loads(response.content)

        # this request should not fail
        if response.status_code != 200:
            logging.critical(
                "connection test failed - got HTTP status \"{}\" with message \"{}\"".format(response.status_code,
                                                                                             content["err_msg"]))
            exit(2)

        logging.info("connection successful - logged in as user \"{}\"".format(content["username"]))

    except Exception as e:
        # if the network connection fails, GalaxyInstance will throw an exception
        logging.critical("exception during connection test - \"{}\"".format(e))
        exit(2)

    # test connection with a GET to /api/whoami

    # TODO do the test
    # resp = gi.make_get_request(galaxy_url+ "/api/whoami")
    # resp = gi.make_get_request(galaxy_url + "/api/configuration?keys=allow_user_deletion").content
    return gi


def get_beacon_histories(gi: GalaxyInstance) -> List[Dict[str, str]]:
    """
    Fetches beacon history IDs from galaxy

        Parameters:
            gi (GalaxyInstance): galaxy instance from which to fetch history IDs

        Returns:
            beacon_histories (List[Dict[str, str]]: List of beacon histories in format {"history_id": "<ID>"}
    """
    beacon_histories_response = gi.make_get_request("{}/api/users/beacon".format(gi.base_url))

    if beacon_histories_response.status_code != 200:
        # TODO error handling
        logging.critical("bla")

    beacon_histories_response_body = json.loads(beacon_histories_response.content)

    if not "beacon_histories" in beacon_histories_response_body:
        # TODO error handling
        logging.critical("bla")

    return beacon_histories_response_body["beacon_histories"]


@dataclass
class GalaxyDataset:
    """
    Representation of a galaxy dataset

    Contains attributes that are used in the scope of this script which is only a subset of
    attributes returned by the galaxy api
    """
    name: str
    id: str
    uuid: str
    extension: str
    metadata_dbkey: str



def get_datasets(gi: GalaxyInstance, history_id: str) -> List[GalaxyDataset]:
    """
    Fetches a given histories datasets from galaxy

        Parameters:
            gi (GalaxyInstance): galaxy instance to be used for the request
            history_id (str): (encoded) ID of the galaxy history

        Returns:
            datasets (List[GalaxyDataset]): list of all datasets in the given history
    """

    datasets = gi.histories.show_matching_datasets(history_id)

    datasets_to_import: List[GalaxyDataset] = []
    for d in datasets:
        if d["deleted"]:
            logging.info("skipping deleted dataset {}".format(d["name"]))
            continue
        if not ["vcf", "vcf_bgzip"].__contains__(d["extension"]):
            logging.info("skipping dataset {} with non-variant format {}".format(d["name"], d["extension"]))
            continue
        datasets_to_import.append(GalaxyDataset(d["name"], d["id"], d["uuid"], d["extension"], d["metadata_dbkey"]))

    # cast datasets from response as Dataset objects
    return datasets_to_import


@dataclass
class BeaconMetadata:
    """
    Dataset metadata to be consumed by beacon during import.

    All metadate fields will be accessible via beacons api. Given values should be ones that can
    shared without exposing too much information (i.e. do not expose internal galaxy IDs, real filenames, etc)

        Attributes:
            name (string): The name for the dataset
            dataset_id (string): Unique identifier for the dataset - any string may be given
            description (string): A short description of the dataset
            assembly_id (string): Identifier of the reference genome used to calculate the variants
            external_url (string): URL will be returned by beacon api along with matching variants from this dataset
            access_type (string): beacon supports CONTROLLED, REGISTERED, PUBLIC
                TODO... other attributes
    """
    name: str
    dataset_id: str
    description: str
    assembly_id: str
    external_url: str
    access_type: str
    create_date_time: str
    update_date_time: str
    call_count: int
    version: str
    variant_count: int

    def __json__(self) -> str:
        """
        Returns metadata json, converting snake_case attribute names to camelCase
        """
        data = {}
        # looping through all attributes and their values
        for attribute, value in vars(self).items():
            # convert attribute name to camel case (replacing _([a-z]) by upper case of the matching letter)
            key = re.sub(r"_([a-z])", lambda x: x.group(1).upper(), attribute)
            data[key] = value
        return json.dumps(data)


def prepare_metadata_file(dataset: GalaxyDataset, output_path: str) -> bool:
    """
    Prepares a metadata file to be consumed by beacon database import scripts

        Parameters:
            dataset (GalaxyDataset): The dataset which will be described by the metadata
            output_path (string): Full destination path for the metadata file

        Returns:
            True if the metadata file has been written and False if ..
    """

    # TODO find assembly ID information
    assembly_id: string

    # Gruppieren nach GRCh38 / 3

    # extract version number from dbkey
    # genome patch level is discarded here
    match = re.match(r"GRCh|hg([0-9]+).*", dataset.metadata_dbkey)
    if match is None:
        logging.warning("found no reference for dataset {} - defaulting to {}")
        #TODO - no reference - no import
        return False

    assembly_id = match.group(1)


    # TODO get call count from dataset
    call_count: int
    call_count = 123

    # TODO get call count from dataset
    variant_count: int
    variant_count = 123

    # assemble metadata from collected information
    metadata = BeaconMetadata(
        name=dataset.uuid,
        dataset_id=dataset.uuid,
        description="some message about where this file comes from",
        assembly_id=assembly_id,
        external_url="usegalaxy.eu",
        access_type="PUBLIC",
        create_date_time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        update_date_time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        call_count=call_count,
        version="v0.4",
        variant_count=variant_count
    )

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(metadata.__json__())


def download_dataset(gi: GalaxyInstance, dataset: GalaxyDataset, filename: str) -> None:
    """
    Downloads a dataset from galaxy to a given path

        Parameters:
            gi (GalaxyInstance): galaxy instance to download from
            dataset (GalaxyDataset): the dataset to download
            filename (str): output filename including complete path

        Returns:
            Nothing

    """
    try:
        gi.datasets.download_dataset(dataset.id, filename, use_default_filename=False)
    except Exception as e:
        # TODO catch exceptions
        logging.critical("something went wrong while downloading file - {}".format(e))


def beacon_import(dataset_file: str, metadata_file: str) -> None:
    """
    Import a dataset to beacon

        Parameters:
            dataset_file (str): full path to the dataset file
            metadata_file (str): full path to a file containing matching metadata for the dataset
                metadata should be in BeaconMetadata format

        Returns:
            Nothing

        Note:
            This function uses BeaconDB from the beacon-python package found at https://github.com/CSCfi/beacon-python

            The connection settings are configured by ENVIRONMENT_VARIABLE or  "default value" if not set

                host: DATABASE_URL / "localhost",
                port: DATABASE_PORT / "5432",
                user: DATABASE_USER / beacon",
                password: DATABASE_PASSWORD / beacon",
                database: DATABASE_NAME / beacondb",

    """
    loop = asyncio.get_event_loop()

    db = BeaconDB()
    loop.run_until_complete(db.connection())
    loop.run_until_complete(
        db.check_tables(["beacon_dataset_table", "beacon_data_table", "beacon_dataset_counts_table"]))

    dataset_vcf: VCF
    dataset_vcf = VCF(dataset_file)

    print(dataset_file)
    variant: Variant
    for variant in dataset_vcf:
        print(variant.num_called)
        print(variant.aaf)
        exit(10)

    # Insert dataset metadata into the database, prior to inserting actual variant data
    dataset_id = loop.run_until_complete(db.load_metadata(dataset_vcf, metadata_file, dataset_file))

    # Insert data into the database
    # TODO !!!! vcf files with no AC info field are not importet - db_load.py:307
    # 1. annotattionsfunktion? 2. file skippen? kommunikation
    #
    loop.run_until_complete(db.load_datafile(dataset_vcf, dataset_file, dataset_id, min_ac=1))


def cleanup(dataset_file: str, metadata_file: str):
    """
    Removes given files
    """
    os.remove(dataset_file)
    os.remove(metadata_file)


def main():
    args = parse_arguments()
    set_up_logging(args.verbosity)
    gi = set_up_galaxy_instance(args.galaxy_url, args.galaxy_key)

    for history in get_beacon_histories(gi):
        for dataset in get_datasets(gi, history["history_id"]):
            # dataset import happens here

            dataset_file = "/tmp/dataset-{}".format(dataset.uuid)
            metadata_file = "/tmp/metadata-{}".format(dataset.uuid)


            # TODO VCF header checken
            # WENN kein AF oder kein AC+AN Feld -> Warnung + Überspringen

            download_dataset(gi, dataset, dataset_file)
            prepare_metadata_file(dataset, metadata_file)

            beacon_import(dataset_file, metadata_file)


if __name__ == '__main__':
    """
    Execute the script
    """
    main()
