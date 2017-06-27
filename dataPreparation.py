import time
import shutil
import os

# Function receiving the input file handle; originally called 'get_fa'
def database_init(variables, db):
    # Initializing redis database

    start_time = time.time()

    pipe = db.pipeline()

    with open(variables.filename) as fIn:

        counter = 0
        for line in fIn:
            # For every header
            if line[0] == ">":
                # read id, template start, template end, query start, query end
                read_id, template_position_start, template_position_end, sequence_start, sequence_end = line[1:].split()
                # query database start and end
                read_position = [int(sequence_start), int(sequence_end)]
                # Beginning and end position within the template
                template_position = [int(template_position_start), int(template_position_end)]

                counter = counter + 1
            else:
                # Sequence is saved (without the last character) in the read dictionary (r)
                read_sequence = line[:-1]

                pipe.hset('read_sequence', read_id, read_sequence)
                pipe.hset('read_position', read_id, 0)
                pipe.hset('template_position', read_id, template_position)

                #pipe.hmset(database_key, database_value)

                # Values are set as a dictionary object and need to be retrieved using the following command:
                # rserv.hget(read_id (e.g. '1100.1'), type_of_entry (e.g. 'read_position))

                counter = counter + 1
            # When counter reaches 1000, execute pipeline
            if counter == 1000:
                # Execute pipeline
                pipe.execute()

                # Start a new pipeline
                pipe = db.pipeline()

                counter = 0

    # Flush the rest of the pipeline
    pipe.execute()

    # Dumping the read database to the disk
    print('Saving original database.')

    db.save()
    # Checking if a directory 'original database' exists and if not creating it
    if not os.path.exists('original_database'):
        os.mkdir('original_database')

    # Copying backup database in the original directory file
    shutil.copy('dump.rdb', 'original_database/dump.rdb')


    elapsed_time = time.time() - start_time

    print('Time it took {} to create a Redis database:'.format(elapsed_time))

    return

def combine_duplicated_reads(variables, db):

    # Loading lua script for database splitting
    with open(variables.LUA_PATH + '/uniqueSeq.lua', 'r') as f:
        uniqueSeqScript = f.read()

    with open(variables.LUA_PATH + '/dereplicate.lua', 'r') as f:
        dereplicate = f.read()

    with open(variables.LUA_PATH + '/addLength.lua', 'r') as f:
        addLength = f.read()

    # Registering lua script to the server and retrieving the pointer
    uniqueSeq = db.register_script(uniqueSeqScript)
    dereplicate = db.register_script(dereplicate)
    addLength = db.register_script(addLength)

    start = time.time()

    uniqueSeq(keys=['read_sequence', 'hashNames', 'dereplicated'])

    dereplicate(keys=['hashNames', 'read_sequence', 'read_position', 'template_position'])

    addLength(keys=['read_sequence'])

    print('It took', time.time() - start, 'seconds to dereplicate sequences.')

    return None

# I don't think this is a necessary database.
# r_pos is used as a legacy position, which is currently saved in read_db_original
# All read_position_db entries are set to 0.

#def initialize_read_pos(variables, dat):
    # This function starts a database with read ids and position 0 for each on of them.
    # Retreaves keys from the read database
#    for read_id in dat.read_db:
        # Sets zero for each of the keys
#        variables.read_position_db[read_id] = 0
#    return None