import time

# Function receiving the input file handle; originally called 'get_fa'
def database_init(variables, read_db, read_db_original):
    # Initializing redis database

    start_time = time.time()

    pipe = read_db.pipeline()

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
                database_value = {'read_sequence' : read_sequence, 'read_position' : read_position, 'template_position' : template_position}
                database_key = read_id

                pipe.hmset(database_key, database_value)

                # Values are set as a dictionary object and need to be retrieved using the following command:
                # rserv.hget(read_id (e.g. '1100.1'), type_of_entry (e.g. 'read_position))

                counter = counter + 1
            # When counter reaches 1000, execute pipeline
            if counter == 1000:
                # Execute pipeline
                pipe.execute()

                # Start a new pipeline
                pipe = read_db.pipeline()

                counter = 0

    # Flush the rest of the pipeline
    pipe.execute()

    # Dumping the read database to the disk
    read_db.bgsave()

    elapsed_time = time.time() - start_time

    print('Time it took {} to create a Redis database:'.format(elapsed_time))


    print('Read database contains {} keys.'.format(read_db.dbsize()))

    print('Read database original contains {} keys.'.format(read_db_original.dbsize()))

    # We are leaving the read_db_original instance as a slave to read_db instance to make sure the transfer of data has
    # time to finish.
    # Since the transfer is happening on the background, we can keep on doing other stuff.

    return

def combine_duplicated_reads(dat):
    # Dereplicates database reads
    sequence_to_read_id = {}
    # Iterating through the read_db and saving seq to a new dictionary
    # If sequence already exists in a dictionary, it's not passed in anymore
    for seq_id, seq in dat.read_db.items():
        if seq not in sequence_to_read_id:
            sequence_to_read_id[seq] = []
        # Read ids of dereplicated sequences are appended to the original key as a list
        sequence_to_read_id[seq].append(seq_id)

    read_db_cleaned = {}
    # Iterating through the dictionary of dereplicated reads
    for seq in sequence_to_read_id:
        # Joining dereplicated reads with the '|' separator
        new_id = "|".join(sequence_to_read_id[seq])
        # Saving in the new database
        read_db_cleaned[new_id] = seq

    dat.read_db = read_db_cleaned

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