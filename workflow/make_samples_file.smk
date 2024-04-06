configfile: "workflow/config.yml"

rule target:
    input:
        config["output_dir"] + "/samples.txt"

checkpoint check_samples:
    input: 
        script = "scripts/check_samples.py",
        raw_data_dir = config["raw_data_dir"]
    params:
        output_dir = config["output_dir"]
    output: 
        config["output_dir"] + "/samples.txt"
    shell:
        f"""
        python {{input.script}} \
        --idat-dir {{input.raw_data_dir}} \
        --output-dir {{params.output_dir}}
        """ 