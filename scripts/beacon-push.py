#!/usr/bin/env python3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORTS
# TODO document dependencies and provide requirements.txt (and make some kind of package)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import asyncio
import asyncpg
import logging
import json

# import utilities from beacon-python
# pip install git+https://github.com/CSCfi/beacon-python
#
#TODO can we install it like this 
from beacon_api.utils.db_load import BeaconDB

# cyvcf2 is a vcf parser
from cyvcf2 import VCF


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOGGING SETUP
# TODO logging should be configurable and support readable and verbose formats
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
formatting = "[%(asctime)s][%(name)s][%(process)d %(processName)s][%(levelname)-8s] (L:%(lineno)s) %(module)s | %(funcName)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=formatting)
LOG = logging.getLogger("beacon-push")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONFIGURATION VARIABLES
# TODO these configuration variables need to be read from command line args and/or os.env
# TODO use argparse
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# base directory from which the galaxy object store can be accessed
galaxy_objects_base_dir = "/home/benno/beacon/docker-galaxy-stable/compose/export/galaxy/database/objects"

# a directory where temporary files can be stored
tmp_dir = "/tmp"

# access to galaxy database
galaxy_database_host = "172.19.0.3"
galaxy_database_port = "5432"
galaxy_database_user = "galaxy"
galaxy_database_pass = "chaopagoosaequuashie"
galaxy_database_name = "galaxy"

# script parameters
dry_run = False


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GalaxyDB class
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TODO use api to fetch objects
class GalaxyDB:
    def __init__(self) -> None:
      self._conn = None
    
    async def connect(self):
        """ establish a connection to the galaxy database """
        try:
          #TODO use a library that can connect to diffrent SQL databases
          #TODO use sqlalchemy
          self._conn = await asyncpg.connect(
            host     = galaxy_database_host,
            port     = galaxy_database_port,
            user     = galaxy_database_user,
            password = galaxy_database_pass,
            database = galaxy_database_name,
          )
        except Exception as e:
            LOG.error(f"failed to connect to galaxy database - {e}")


    #TODO get files via api (use bioblend + download file)
    async def getDatasetUUIDs(self) -> list:
        """ Fetches all vcf and vcf.gz files that need to be pushed to beacon. Returns a list of their UUIDs."""
        #TODO can there be a better filter instead of name+tag combination?
        #TODO is it worthwhile to supply tag ID and save a join?
        #TODO join user and use beacon setting
        records = await self._conn.fetch(
            """
            SELECT
              d.uuid as uuid
            FROM history h
            INNER JOIN
              history_tag_association ta ON h.id = ta.history_id
            INNER JOIN
              tag t ON ta.tag_id = t.id and t.name = 'beacon-share'
            INNER JOIN
              history_dataset_association da ON h.id = da.history_id AND da.extension IN ('vcf', 'vcf_bgzip')
            INNER JOIN 
              dataset d on d.id = da.dataset_id
            WHERE
              h.name = 'Beacon Shared Variants'
            """
        )
        return [r['uuid'] for r in records]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions used in main loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def path_from_uuid(uuid: str) -> str:
    """ returns the absolute path of the file with the given UUID """
    #TODO can the path always be generated like this or is this a special case?
    return f"{galaxy_objects_base_dir}/{uuid[0]}/{uuid[1]}/{uuid[2]}/dataset_{uuid[0:8]}-{uuid[8:12]}-{uuid[12:16]}-{uuid[16:20]}-{uuid[20:]}.dat"

def create_metadata_file(uuid: str) -> None:
    """ creates a metadata file for the given uuid """

    #TODO retrieve callCount and variantCount somewhere (galaxy db? file itself?)
    #TODO make a nice description about where to get more information abount the variant
    #TODO should the internal sample uuid kept secret?
    #TODO createTime and updateTime - what can we put here?
    
    with open(f"{tmp_dir}/metadata-{uuid}.json", "w") as metadata_json_file:
        data = {
            "name": uuid,
            "datasetId": uuid,
            "description": "some message about where this file comes from",
            "assemblyId": "i have to find out what this is",
            "externalUrl": "usegalaxy.eu",
            "accessType": "PUBLIC",
            "createDateTime": "2013-05-02 12:00:00",
            "updateDateTime": "2013-05-02 12:00:00",
            "callCount": 123,
            "variantCount": 456
        }
        json.dump(data, metadata_json_file)
        #TODO tmp dir cleanup


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main program
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
async def main():
  #TODO parameter handling should go somewhere here
  gdb = GalaxyDB()
  await gdb.connect()
  datasetUUIDs = await gdb.getDatasetUUIDs()

  db = BeaconDB()
  await db.connection()
  await db.check_tables(["beacon_dataset_table", "beacon_data_table", "beacon_dataset_counts_table"])


  # main loop
  for uuid in datasetUUIDs:
      datafile  = path_from_uuid(uuid)
      metadata_file = f"{tmp_dir}/metadata-{uuid}.json"

      if (dry_run):
          print(datafile)
          continue

      create_metadata_file(uuid)

      vcf = VCF(datafile)
      
        
      # Insert dataset metadata into the database, prior to inserting actual variant data
      dataset_id = await db.load_metadata(vcf, metadata_file, datafile)

      # Insert data into the database
      #TODO !!!! vcf files with no AC info field are not importet - db_load.py:307 
      # 1. annotattionsfunktion? 2. file skippen? kommunikation
      await db.load_datafile(vcf, datafile, dataset_id, min_ac=1)


  # Close the database connection
  await db.close()


#TODO ask about asyncio - why are they using it?
asyncio.run(main())