pipeline {
  agent any
  stages {
    stage('First stage') {
      steps {
        sh 'printenv'
      }
    }
    stage('HD832') {
      steps {
        sh """
              source /etc/bashrc
              virtualenv venv -p python3.8
              source venv/bin/activate
              pip install -r requirements.txt
              mkdir -p HD832
              cp .tests/release/units.tsv HD832/units.tsv
              cp .tests/release/samples.tsv HD832/samples.tsv
              cp .tests/release/resources.yaml HD832/resources.yaml
              cp -r config HD832/

              module load singularity
              module load slurm-drmaa
              snakemake -s workflow/Snakefile --profile .tests/release/profiles -d HD832
           """
       }
    }
  }
}
