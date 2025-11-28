# Changelog

## [1.0.0](https://github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.24.0...v1.0.0) (2025-11-28)


### ⚠ BREAKING CHANGES

* use config and profile to start ctDNA instead of samples.tsv

### feat\

* use config and profile to start ctDNA instead of samples.tsv ([89ed27c](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/89ed27c139c48713530ddde585f0b1877457be0d))


### Features

* add CTAT-splicing ([b95dc9e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/b95dc9e557cd8c98b9e1fcd74bdbce7511e906ca))
* added ctat-splicing ([db4f3e4](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/db4f3e48878ec733c011c753c6b87e466c70f333))
* added miarka build and ref files ([b9a2e6e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/b9a2e6eeb3b500a906dc7c8ccdfa1d67429f5b03))
* added updated CTAT references files ([d2003f9](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/d2003f917ef1883991680654fd2b2393abf0105d))
* handle tc in percent instead of fractions ([d1caa63](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/d1caa63653d4f67962c872e4d79f4e9bd0668422))
* handle tc in percent instead of fractions ([c9b7148](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/c9b7148931d654f2adb2edd9516f60de108425e9))
* use new fusion module version with ctat ([80a265e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/80a265e653118f9ae83448f0da1c1bc4bf77a1ff))


### Bug Fixes

* **build_conda.sh:** Set pipeline version in path ([05b0f06](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/05b0f06de601dee934edc8c12190d441399d93c9))
* **build:** Optional to re-build ref. and apptain. ([23192a5](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/23192a5686bae10d430a1bc8360197410c9176bb))
* correct output file ([16d8759](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/16d8759b7c50873192baad44c9752617f7ffc602))
* correct output file ([7631485](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/7631485d091a3b3c374ffecbdfe2ed226eaccd8d))
* make prio analysis work dynamically ([d8dac01](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/d8dac0154d66939b7b12601892598367a89fc5bc))
* make prio analysis work dynamically ([d53441d](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/d53441d9848636fb7cfc8d6abe81aa79763a3d9f))
* Miarka purecn ([aa25b54](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/aa25b5462cf5eb4f9515214676abb8102f860c67))
* remove old files ([bb2f7e0](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/bb2f7e00ca1685fcb4201ced5255806c09fb73e1))
* **resources.yaml:** Resources optitype ([92dba5a](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/92dba5affb83c97b10b907af8849f029c38e21f2))
* **resources.yaml:** Resources optitype ([3cdd3e2](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/3cdd3e21c34b69467e6bfaa6e31a55b81c967f44))
* rm hla optiype analysis ([13d5b82](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/13d5b82ac224ab2f2611001685b476ef00bed2d7))
* rm hla optiype analysis ([b934e16](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/b934e169d9b3844a11bb9eea39a98b95718dc589))
* rm large unfiltered ITD results ([7689af4](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/7689af4feec7ab761fe03738cfa1d4f679527e60))
* rm large unfiltered ITD results ([cbe9143](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/cbe914351e9f5478db4dc3145452e20d531b3fda))
* Update build_conda.sh ([426f182](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/426f18221d09856d313df2dd06d2587955d8cbd5))
* Update build_conda.sh ([654fefb](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/654fefbbd9043ce4469cc73e43225ed877a7e797))
* Update build_conda.sh ([296e01c](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/296e01cf4f68d98b337327b204f6bcb1548b89a8))
* Update build_conda.sh ([3ec795b](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/3ec795bd41d13607889173b618a524ae048a37b0))
* Update build_conda.sh ([eb51270](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/eb5127017c9fb0172c3688f6fcef0ca87f097976))
* Update design_files.hg19.yaml ([30e663b](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/30e663b4f421818082a0e8806a111af087b2e59d))
* Update design_files.hg19.yaml ([eadb33e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/eadb33e78595445caa556843664707678dc994a3))
* Update resources.yaml ([da925df](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/da925dfe7dfa213b910f5e5776ca6f84c3bd5cc2))
* Update singularity.schema.yaml ([ee6ae08](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/ee6ae088a676cf9a7db887cd5493e8eed036c4ff))


### Documentation

* **README:** Updated build instructions. ([86e66f0](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/86e66f04e9f322f24a9f9b22ef689d162ffd7935))
* Update README.md ([97aeb16](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/97aeb169ea8d782141fb842b487a435a4059242b))

## [0.24.0](https://github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.23.0...v0.24.0) (2025-10-17)


### Features

* cov_and_mut based on umi when umi is used ([bf59619](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/bf59619d0831b136d7fd24332afddaf69fec8e6b))
* cov_and_mut file with aa one letter code ([0b15bb1](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/0b15bb123aec209211e6f0c8ffbde8e8119a527b))
* cov_and_mut file with aa one letter code ([b88f617](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/b88f6173d30bd46c9212cdf74a7700c6f3bfa5ad))
* Coverage_and_mutations files for haloplex ([568ab9e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/568ab9e264b9a90b1726ef11b5c324011ed097b2))
* updated coverage and mutations files ([3dda8c9](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/3dda8c975c352f850cd0e02fd5bf32700a999161))
* updated coverage and mutations files ([c392599](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/c392599a8a388ab655bc42fbfa7e2efe33891e2b))


### Bug Fixes

* adjust GRCh37 VEP cache download location ([#667](https://github.com/genomic-medicine-sweden/Twist_Solid/issues/667)) ([a1eaecb](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/a1eaecb9335edb774d0be4d04203e8997c088b6e))
* fix broken annotation in Mutation_Lung ([09d01c4](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/09d01c44f44bfaeb4f85d247be68dc17cac795e1))
* fix broken annotation in Mutation_Lung ([9f43a82](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/9f43a82c6d58175b5181791b772f41f19d4ba028))
* Gist cov and mut file ([3794e29](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/3794e29c38c899909528916bae93e070a6ddcd8a))
* output file variable name ([164050e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/164050e7af6bf5fe369191929053a23c48a10f3f))
* reference files ([7c1c2d8](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/7c1c2d8a87e42b11a6b6d7550f2b10c924f5691b))
* reference files ([8aea148](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/8aea1480a836bccee24571f48388032b38bf0794))
* Update general_report.yaml ([d295f16](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/d295f16efea18a5e44118c4789fd40157136f93d))
* Update general_report.yaml ([7ce88b9](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/7ce88b944bb658f1b4dccbf7fd171279ef8c6089))
* update mod cov and mut files ([3db8110](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/3db8110063941fa93255c218775ca1030fedb6c9))
* update mod cov and mut files ([1fdc09e](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/1fdc09e441cee1ab7566ec0c5a6452ec5317fb4e))
* Update output_files.yaml ([bde39b8](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/bde39b897297c291a758af52a4f19de3e647666c))
* Update output_files.yaml ([aab581c](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/aab581ca27a9839bba9cdfc649c4c437c28b6ee6))
* Update Snakefile ([8db0fd0](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/8db0fd024d06f1d5e58ab205ae558224d99bc5bb))


### Documentation

* added clarifying info ([e997a71](https://github.com/genomic-medicine-sweden/Twist_Solid/commit/e997a715284255a2607ccd74807880b9060db05c))
