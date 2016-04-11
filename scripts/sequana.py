import os
import shutil
import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """%s --init fix_removal\n""" % prog
        usage += """Examples:

            sequana --init <sequana pipeline>
            sequana run +snakemake options
            sequana report

        """
        super(Options, self).__init__(usage=usage, prog=prog)
        self.add_argument("--init", dest='snakefile', type=str,
                required=True, 
                help="Get the snakefile and config file and possible other files")


def main(args=None):
    from sequana.snaketools import modules, Module, SequanaConfig
    import shutil
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        print(user_options)
    else:
       options = user_options.parse_args(args[1:])

    try:
        module = Module(options.snakefile)
    except:
        print("Invalid module name provided (%s). " % options.snakefile) 
        print("Here are the valid modules from sequana: \n - %s" %
"\n - ".join(sorted(modules.keys())))

    # the module exists, let us now copy the relevant files that is
    # the Snakefile, the config and readme:
    shutil.copy(module.snakefile, "Snakefile")
    shutil.copy(module.config, ".")
    shutil.copy(module.readme, ".")

    # we also crate the directory to hold the raw data
    try:
        os.mkdir("fastq_raw")
    except:
        pass

    # a running script

    with open("sequana.sh", "w") as fh:
        fh.write("#!/usr/bin sh\n")
        fh.write("snakemake -s Snakefile --stats stats.txt -p")
    #os.chmod("sequana.sh", 744)

    print("""You can run snakemake yourself or type::
    sh sequana.sh
    """)


    # figure out from the config file if any files are required
    try:
        # stop here if the file does not exists
        configfile = os.path.split(module.config)[1]
        config = SequanaConfig(configfile)
        if 'requirements' in config.parameters.keys():
            from sequana import sequana_data
            for filename in config.parameters.requirements:
                try:
                    print('Getting %s ' % filename)
                    shutil.copy(sequana_data(filename, "data"), ".")
                except:
                    print("This file %s will be needed" % filename)

    except Exception as err:
        print(err)
        print("No config file found")



if __name__ == "__main__":
   import sys
   main(sys.argv)

