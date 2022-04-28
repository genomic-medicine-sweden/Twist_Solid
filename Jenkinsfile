pipeline {
  agent { label 'marvin-clone' }
  stages {
    stage('HD832') {
      steps {
        sshagent(['jenkins']) {
          sh """
                source /etc/bashrc
                virtualenv venv -p python3.8
                source venv/bin/activate
                pip install -r requirements.txt
                mkdir -p HD832/slurm_out
                cp -r config HD832/
                cp .tests/release/units.tsv HD832/units.tsv
                cp .tests/release/samples.tsv HD832/samples.tsv
                cp .tests/release/resources.yaml HD832/config/resources.yaml
                cp .tests/release/config.yaml HD832/config/config.yaml

                module load singularity
                module load slurm-drmaa
                snakemake -s workflow/Snakefile --profile .tests/release/profile -d HD832 -w 60
             """
         }
       }
    }
  }
}
