from reagoClasses import JSONObject
import os

def parseInput(args):
    # Retrieving an input file

    # Creating a json object containing user settings
    variables = JSONObject()

    try:
        variables.read('settings.json')

    except:
        print("No settings file was found.")

    ##############################################################################################
    # Input file
    # Default value of the 'input' argument is 'None', given that user didnot specify input, let's check settings file.
    if args.IN == None:
        try:
            # If settings file contains input file, we're good to go.
            variables.filename
        except:
                print('ERROR: Please specify an "input file" path using token -i/--input.')
                quit()
    # User have specified input file.
    else:
        variables.update({'filename':args.IN})

    # Checking whether input file exists
    if os.path.isfile(variables.filename) == False:
        print('ERROR: Input file', '\"' + variables.filename + '\"','does not exist. Please double check the path and filename in the -i/--input token or the \"settings.json\" file.')
        quit()

    ##############################################################################################
    # Output folder
    # Default value of the 'output' argument is 'None', given that user didnot specify input, let's check settings file.
    if args.OUT == None:
        try:
            # If settings file contains output file, we're good to go.
            variables.output_dir
            # Adding a slash behind the output folder should it not already contain one
            if variables.output_dir[-1] != "/":
                variables.output_dir += "/"
        except:
                print('ERROR: Please specify an "output directory" path using token -i/--input.')
                quit()
    # User have specified output file.
    else:
        variables.update({'output_dir':args.OUT})
        # Check if output contains slash and if not add it
        if variables.output_dir[-1] != "/":
            variables.output_dir += "/"

    # Checking whether input file exists and if not, creating it.
    if os.path.isdir(variables.output_dir) == False:
        os.mkdir(variables.output_dir)

    # Formating the output
    variables.update({'graph_path':variables.output_dir + 'graph.data'})
    variables.update({'plot_dir':variables.output_dir + 'plot/'})
    variables.update({'rj_dir':variables.output_dir + 'rj/'})
    variables.update({'full_genes_path':variables.output_dir + 'full_genes.fasta'})
    variables.update({'fragments_path':variables.output_dir + 'fragments.fasta'})

    # If readjoiner (rj) directory doesn't exist yet, it will be generated.
    if os.path.exists(variables.rj_dir) == False:
        os.mkdir(variables.rj_dir)

    ##############################################################################################
    # Retrieving user input values and renaming them.
    if args.READ_LENGTH == None:
        try:
            # If settings file contains read length, we're good to go.
            variables.READ_LEN
        except:
                print('ERROR: Please specify a "read length" using token -l/--read_length')
                quit()
    else:
        variables.update({'READ_LEN':int(args.READ_LENGTH)})

    ##############################################################################################

    if args.OVERLAP == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.MIN_OVERLAP
        except:
                print('ERROR: Please specify an "overlap" using token -ol/--overlap')
                quit()
    else:
        variables.update({'MIN_OVERLAP':variables.READ_LEN * args.OVERLAP})

    ##############################################################################################

    if args.TIP_SIZE == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.TIP_SIZE
        except:
            print('ERROR: Please specify a "tip size" using token -t/--tip_size')
            quit()
    else:
        variables.update({'TIP_SIZE': int(args.TIP_SIZE)})

    ##############################################################################################

    if args.CONFIDENCE_BASE == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.CONFIDENCE_BASE
        except:
            print('ERROR: Please specify a "confidence base" using token -c/--confidence_base')
            quit()
    else:
        variables.update({'CONFIDENCE_BASE': int(args.CONFIDENCE_BASE)})

    ##############################################################################################

    if args.ERROR_CORRECTION_THRESHOLD == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.ERROR_CORRECTION_THRESHOLD
        except:
            print('ERROR: Please specify an "error correction threshold" using token -e/--error')
            quit()
    else:
        variables.update({'ERROR_CORRECTION_THRESHOLD': int(args.ERROR_CORRECTION_THRESHOLD)})

    ##############################################################################################

    if args.FULL_LENGTH == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.FULL_LENGTH
        except:
            print('ERROR: Please specify an expected "full length" using token -f/--full_length')
            quit()
    else:
        variables.update({'FULL_LENGTH': int(args.FULL_LENGTH)})

    ##############################################################################################
    # todo : Figure out whether this parameter is actually used anywhere.

    if args.PATH_FINDING_PARAMETER == None:
        try:
            # If settings file contains minimum overlap, we're good to go.
            variables.PATH_FINDING_PARAMETER
        except:
            print('ERROR: Please specify a "path finding parameter" using token -p/--path_finding')
            quit()
    else:
        variables.update({'PATH_FINDING_PARAMETER': int(args.PATH_FINDING_PARAMETER)})

    ##############################################################################################

    # todo : Figure out what that is
    variables.NEED_DEFLANK = False

    variables.LUA_PATH = 'lua_scripts'

    #dat.update({'read_db':{}})

    #variables.update({'r_pos':{}})
    #variables.update({'cm_pos':{}})
    #variables.update({'read_db_original':{}})
    #variables.update({'read_position_db':{}})

    variables.write('settings.json')

    #import pdb; pdb.set_trace()

    #TODO g.MIN_OVERLAP = args.READ_LENGTH*args.OVERLAP
    #TODO g.TIP_SIZE = args.TIP_SIZE
    #TODO g.CONFIDENCE_BASE = args.CONFIDENCE_BASE
    #TODO g.ERROR_CORRECTION_THRESHOLD = args.ERROR_CORRECTION_THRESHOLD
    #TODO g.READ_LEN = args.READ_LENGTH
    #TODO g.FULL_LENGTH = args.FULL_LENGTH
    #TODO g.PATH_FINDING_PARAMETER = args.PATH_FINDING_PARAMETER
    #TODO g.NEED_DEFLANK = False

    #TODO g.graph_path = g.output_dir + "graph.data"
    #TODO g.plot_dir = g.output_dir + "plot/"
    #TODO g.rj_dir = g.output_dir + "rj/"
    #TODO g.full_genes_path = g.output_dir + "full_genes.fasta"
    #TODO g.fragments_path = g.output_dir + "fragments.fasta"

    return variables