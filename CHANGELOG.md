# Changelog

## [0.7.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.6.1...v0.7.0) (2023-08-17)


### Features

* add result file for combined purecn and pathology ([fb17cc2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fb17cc25f112994a8a7047f56072d0a0ae32ac75))
* added duplication % to multiQC ([2923ebf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2923ebf492b672993d71afa992f9d577b0879149))
* added picard mark duplicates of bam-files for QC ([0f33656](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0f33656a54c007d027591d3e9cf117ca0637b148))
* added read group function for STAR ([c224c17](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c224c177eeb72d8fbbc3d42c7af1cc1ec3068787))
* added RG to STAR and changed bam file for QC ([1f68585](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1f685859ae5e36716a5b3e5716de3f90dce0c4f7))
* added rule for modifying MBQ in vcf ([72dc309](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/72dc309cebb3eacb6d68f1c13fe65e3f904d05f2))
* added rule for modifying MBQ in vcf ([6c26626](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6c26626f45857a00eea1e8f257cd9370ee1163c9))
* added rule for modifying MBQ in vcf ([188f418](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/188f418a5385fad636e7603a3baf14eda4ee41c2))
* change pureCN cutoff to 0.35 ([4e06643](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4e06643375303ff8718607e42da314fe1bed6356))
* choose purecn if tc > 30% and pathology otherwise ([20292ab](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/20292abf636c1514cb4c3d5e8a6c395bc36c0ef4))
* harder filtering ([7cae319](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7cae319b8ed460a194788ecf2f81eb0441c010ea))
* make two tsv reports using different gene lists ([0737c60](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0737c605dea5a90d87735cc0dd2b97d0f32882f5))
* test_input_all.tsv for v0.7.0 ([27df350](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/27df350071e9d23f4a024544d79f82fe6498d16a))
* test_input_VAL2022.tsv for v0.7.0 ([26f5157](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/26f515742d0d1587b57acf0cda4ba2dd0d526a88))
* use filtered vcf with both germline and somatic variants ([f6c5cc3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f6c5cc3ac86c85e40130edf2edfac37da1379966))
* use gatk2 for purecn ([ee2e2bc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ee2e2bc98e034cfc01f405adc9e35cf09ff03e55))
* use germline vcf for purecn ([eb0eaf5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/eb0eaf5c2bbd26c5e3f960488e147c65092ae02a))
* use purity file directly from purecn to also get ploidity ([b397a3c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b397a3c9a486d8ed1d8b0f326ae553d088c3cb9d))
* use vaf and snv filtered vcf with both germline and somatic variants ([48f3563](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/48f35631991ede274ee5daeea30a9ba252964333))


### Bug Fixes

* add germline flag to vcf ([1e8de1b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1e8de1bb1c4741b61ab4915298f82dd3c0ed4bdc))
* add missing filter tag ([7c01e2d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7c01e2d17205b03938ea73d6afbdbf412fea3238))
* annotate using missing sites instead ([de172b0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/de172b06028c2acab5b0f82a9ce9d699cb35c9f8))
* bug fixes ([224a6a4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/224a6a48a9f3db711a9c86f4789235619557fce4))
* change checkpoint to rule ([5e04c9a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5e04c9acad4c27049b14ae1829d4eb840fa0fc5f))
* change path to new normals ([7a3b420](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7a3b4209d7d16ba8672e464f0845c1b4a4ad6ae6))
* correct header in cnv report file ([c947152](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c947152845a43019db577035d549f9b7fe55c5d5))
* correct output name for purecn reference ([92af141](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/92af1410b325942fd0a12c04c7b971917c2c758a))
* correct rule import from wrong module ([275a60d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/275a60d0ed2fe7c8a86f92d0882c4dfd29fd90a4))
* delegate schema validation to reports module ([b5dada0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b5dada0172b5b94e7eb89c57e5b164ef0fae29d8))
* do not filter large cnvs based on frequency in database ([4b86637](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4b86637cf673856a7b1899c0c2b27db635fcec46))
* get correct tc to html report ([c394cba](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c394cbabefd5f4443c37f53bf21473ed4a6a67a9))
* handle empty purecn file ([df89fe5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/df89fe570f1651ef8a6a08180309144e71d51e22))
* import spelling mistake ([f3e7248](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f3e72484c31130921ff5c860e836139c266cd024))
* moved result file to additional files ([c6c7574](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c6c75742c39138a16c491c2795134f2982bf7fbd))
* properly overrule the `get_tc` function ([701ceb8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/701ceb8f3a381e71c43874dc417cba9838a24b0f))
* purecn_modify_vcf bugfix ([e910e18](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e910e18b8588e2204622b21536b4595160fb7846))
* redefine rule to use new params in config ([5ef5dac](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5ef5dac3536a90d26320330058c98c9ba81b9f38))
* return correct tc ([45485c7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/45485c76c45bdd96ea87f58c39904f30a608ef1c))
* solve different wildcards in rule error ([62082a1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/62082a1a004676e4b0b38ba9b7d09ef6c18e204b))
* spelling error of Exception ([8d6fbf9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8d6fbf9c8916b79987b598242e8ebe3ec2ab681d))
* tabix of annotation database ([1f70635](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1f70635d33a2c05d2d32b2c07e2f82286c3ae817))
* use correct genome ([66e3c8b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/66e3c8b04f2e0457117d2c9493e86eff071c6c92))
* use correct get_tc ([ba5272e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ba5272e5c86ec743c1984c1834b5b710d3ea02aa))
* use correct interval file ([92c0e41](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/92c0e419069c174c817eea1161088b288240558b))
* use Illumina for platform ([130f3a2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/130f3a2d8160ef723c3d184c3981b132bbc30764))


### Documentation

* add readthedocs link to readme ([7c61ef3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7c61ef3d2c8ac8dd641ac4a01fc4bbcfbcc979a8))
* update CNV HTML report documentation ([3478612](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/34786124f7cb071868708268b7a99dd35da09228))
* update mention of output spec in docs ([60b9800](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/60b9800cf3b6549b1fd53cebd08afee2c643f2c8))
* update readme ([8673794](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/86737947780e45af1d4fa0e984ad4e142d43950c))
* updated cnv documentation ([a503ddd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a503dddfd4e49e16401e0d25029de0fc4654c4f8))
* updated cnv documentation ([ef6786a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ef6786ab7bf09119381ca27bf2de88c1c46bb6b0))

### [0.6.1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.6.0...v0.6.1) (2023-04-28)


### Bug Fixes

* fix building of read the docs ([8043566](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8043566e569df32fb6a9859f0a61bc3c161f95f0))


### Documentation

* moved build yaml ([a7a3bfe](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a7a3bfe1a3917786c0ef72e86d84ee6d3ab1daf3))
* update build yaml ([c3455aa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c3455aa5bda22a47d14c49634807de4403665a00))
* update build yaml ([aad3220](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aad32201832060c1dd9d2ea936d8aa13c147b05e))
* update build yaml ([6173387](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6173387939d6897243e300017456799f46163c67))
* update build yaml ([4830711](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4830711e23784f4a724ad95e827e9c8d7f645aff))

## [0.6.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.5.0...v0.6.0) (2023-04-28)

## Features

###Documentation
First relase of read the docs (https://twist_solid.readthedocs.io/en/v0.6.0/)

### CNV
- New small cnv amplification caller (in house script) for the 20 relevant amplification genes which is very similar to the small deletion caller. 
- The CNV.tsv report now also includes amplifications called by the small amplification caller

## Bugfixes
- New lines added to variants in the TMB report

## Changes in config.yaml
- Small amplification caller: new in config
- CNV tsv report: added option for small amplification filtering

## Hydra modules with releases
- prealignment: v1.0.0 (No change)
- alignment: v0.3.1 (No change)
- snv_indels: v0.3.0 (No change)
- annotation: v0.3.0 (No change)
- filtering: v0.1.0 (No change)
- qc: v0.3.0 (No change)
- biomarker: v0.3.1 (fixed TMB report bug)
- cnv_sv: v0.3.1 (No change)


### Features

* added caller for small amplifications ([61181d5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/61181d51462458c850813641e13c6d9aba2685fe))


### Bug Fixes

* add missing parameter to call of create_tsv_report ([36f3508](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/36f35080259d5edd1f48cb7afc33d72bc0964fd3))
* bugfixes and config tweaking ([e8313ff](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e8313ff7f9165d5838d791247ab399ca3bfb9dda))
* match rule input names ([981d2a6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/981d2a6525cde093754a8f4c6e1808f2bd15e153))
* parameter tweaking ([16ab328](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/16ab3284af92067c4fcbc618c6a9a1f4eebc1294))
* update fusions and biomarker tags ([9dfadef](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9dfadef44e4f7da4397095bd1051cff639162e61))


## [0.5.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.4.0...v0.5.0) (2023-04-19)
# Release notes
For more details on features and bug fixes see further down.

## Features

### CNV
- The CNV.hmtl report now reports TC content
- The CNV.hmtl report now reports VAF values in the variant table
- The CNV.tsv report now also includes deletions called by the small deletions caller

### TMB
- Improved TMB-calculations. Finds more true variants and have better correlation compared to TSO500. Does not use any panel of normals anymore making the calculations more independent on sequencing platform.

### DNA fusions
- Added DNA fusion calling using Fuseq-WES with superior results compared to GeneFuse
- GeneFuse: Added filtering of the ERG gene

### RNA exon skipping
- Only report MET exon 14 skipping and EGFRvIII and not other potential skipping events in these genes

## Bugfixes
- Copy .bai file with timestamp instead of creating it so that it is not removed by snakemake

## Changes in config.yaml
- TMB: new and updated config options for tmb rule
- FuseqWES: Added config for fuseq_wes rule
- FuseqWES filtering: Added config for filter_fuseq_wes rule

## Hydra modules with releases
- prealignment: v1.0.0 (No change)
- alignment: v0.3.1 (No change)
- snv_indels: v0.3.0 (No change)
- annotation: v0.3.0 (No change)
- filtering: v0.1.0 (No change)
- qc: v0.3.0 (No change)
- biomarker: v0.3.0 (TMB updated with more config options)
- cnv_sv: v0.3.1 (No change)

### Features

* add BAF to CNV table ([8f08da1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8f08da1e2bfad41b838c21778f489305117ce475))
* add tumor cell content to report ([3740ab6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3740ab6eee7d647d45b75728b05297af2c677efa))
* added fuseq_wes rule and filtering to pipeline ([3702411](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/37024115e0786323696858da45cd954c9f45a466))
* added small cnv deletions to tsv report ([2bd67dd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2bd67dda576d98ca5b73f48cecbf08c248ac9df2))
* added trancript black list file ([f4d3aaa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f4d3aaa07a7d4fa6b2b7252d2b90f2187a17b7bd))
* change purity result file to purecn output for additional info ([270718e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/270718eef1d9f7f9672e89ccfcd24b80a32682f6))
* ERG filtering ([9d57325](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9d57325d550a3f2abeeae593e483d1f1fbc69228))
* improved tmb calculations independent of artifacts and background ([452ca47](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/452ca478972576704e9d2e280932f7dde271dc40))
* move non-essential result files to additional files ([7883b8e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7883b8e38d82ef360b02dc2a9c0b1a612febbcf9))
* move small deletions to additional result files ([bd71f95](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/bd71f955e67ceca96559b5793a728c5016e5df16))
* report only MET exon 14 and EGFRvIII ([e6865d0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e6865d02364aca8bd5db7c97ea05395db7d0998e))
* **script:** add flag for noisy fusions ([ed5b5d8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ed5b5d880180794d7fba7c33d93c9d464bef56e5))
* update to latest hydra fusion develop tag ([9bcd34b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9bcd34ba4e6c9b3b4e858766a184026d998a1f3c))
* update to new reference module tag ([3844f52](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3844f52240d2e44c395465f84d631912098bebef))
* updated to biomarker develop tag ([860ebf6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/860ebf6049d076d90a06fca7020d619d84bd3524))


### Bug Fixes

* better table column names ([6419696](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6419696d245af09b4a1221fa6baf70cc4f102c41))
* bugfix ([284520e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/284520eee3ebb9b1fc874df254b24db037e0898f))
* **config:** purecn rules ([aee6424](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aee642403369ae375f9b5f4b07595e262867de76))
* **config:** variable name changed in new tag ([d01fc93](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d01fc93bf78ce9ad37ddaf9383b8fbbbd83d4abf))
* copy bai files instead of creating ([feecb70](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/feecb705fe30e1d9223ffe4206156a73fe82b606))
* docker version ([198f648](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/198f6481acc8c2b4528965b90f741be920965d5c))
* docker version ([f022c57](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f022c57c8169c4fcac2b207e1e701b76987d8d63))
* fix incorrect annoation order of BRAF ([5f65a3d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5f65a3d4246007e8c67a7626ca8bd1023db1f1cd))
* input file name fix ([f1c0694](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f1c06941a64cfce877af4f2e8ade28b5976ee481))
* preserve timestamps when copying ([ff1b17e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ff1b17e8d5abc10db38f7e00d905588b7c9f5a88))
* ruleorder ([21ba998](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/21ba99867428bdbdb491bfbc93545412124327ea))
* **schema:** update variable name ([c3ffb2d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c3ffb2d1c26b5c7eb5836e5dfe7c84c18d58e86d))
* tag name ([1fa43ab](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1fa43abe9f627b65ac6ef8c71ba0805bb7609700))
* update biomarker tag ([518e738](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/518e738d66b750ff2f2f2b31ab9341a1d6b5ddb3))
* Update config/config.yaml ([daa344e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/daa344e82e49b7a87eff022abedbd1a31b87a419))
* vaf and annotations plotted out-of-bounds ([44b8ce8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/44b8ce8b57241e1c0056cd967de8490080f17199))


### Performance Improvements

* increase allowed memory usage for fuseq wes ([10f9804](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/10f98047510248b730408c97834fe5f23f30f96b))


### Documentation

* add missing tumor_content column to sample schema ([7dbfa63](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7dbfa631b46b55a46e731630d5ead6acaa566b5e))

## [0.4.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.3.0...v0.4.0) (2023-02-07)


### Features

* adapt rules and report to new structure ([6990d24](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6990d24562db687e1aaa16280f33e5a57c5d9851))
* add gridlines to genome view ([b9df3b5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b9df3b5a21fcedd5b57506b492548bc2ed1343e8))
* add guides to chromosome view ([06b8cbf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/06b8cbfe0f2b0c0730cd3c1a2bcc166127f04f51))
* add metadata to report ([198644c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/198644c3f330970bba25d370600f4efa426ca4fc))
* add script for merging JSON files ([4008cc8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4008cc8e353104a7d092e04ed01f47e54f7960aa))
* added arguments and files from Erik ([f52961c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f52961cf97ed8f6f1e48a8fa44563d7d247a4df8))
* added blacklist for bad VAF regions ([27f4201](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/27f4201ffbcb20645942bf32065aa828c82eb82e))
* added blacklist to cnvkit_call and cnvkit_scatter ([3003011](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/300301175cbf0838c8cced34648814cf99fa0341))
* added header info in vcf warning of the AD change for QCI ([771322e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/771322ee604c2fca2e55268c5eaadcbc47cb2167))
* added purecn reference generation ([fe4ac35](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fe4ac359abe21ed54f36f5b14f6c0662bf297b0f))
* added purecn tc estimation to results ([4214001](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/421400124c62134101b5e5feb200c215433ff171))
* added purecn to pipeline ([80d3bf5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/80d3bf561db5f7cc2308e84aef0b352329261905))
* Added rule call_small_cnv_deletions calling smaller deletions in cnv data from GATK ([0ac2f58](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0ac2f589d8e62971bbcb869a2135127d1cdeb8a0))
* added rule to fix vcf AD values that is used in QCI to calculate AF ([edb9769](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/edb976978a1af5d7fb6ccaa8751bdf2b78a9fd00))
* added total supporting reads to fusion report ([0f79608](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0f796084e91887dde1a01c115770135854a5154c))
* changed to internal segmentation ([00661ca](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/00661ca01492be079dcaaeac7c8ddda1cc031318))
* convert output from single caller to JSON ([6d70961](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6d70961f16ffcc33b93fed74aa3e1c29c91c38a2))
* create separate print layout ([86da7d1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/86da7d1040a67bd41a453a3bc9a60cc2674f289d))
* filtering of loh cnvs also using baf to keep copy neutral cnvs ([d55e470](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d55e470db7138390e2d1f097e0ce5431a645e09a))
* import function through the module ([8af029d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8af029d2a5034627c929866203a248678e4914d9))
* indicate points outside zoom range ([24745d7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/24745d7339a1d5107369d8a22a5e683bbb5bfa32))
* more flexible layout ([9cfb79d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9cfb79d26c225d861e95efe8f27314cb2de8f9f3))
* remove high VAFs from cnvkit plots ([558e55b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/558e55b0a3d212a8d75cdd0a5575ecff54b1c4c1))
* run purecn with GATK cnv ([6470af4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6470af4cc219fbff9bdb2c178c851499456a51db))
* scale y-axis to extent of data when zooming ([39f81d3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/39f81d3fb97fd1d04636daae4c5257d584dc0e5d))
* show results from other callers on hover ([73f7f88](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/73f7f88e78ca3d61f504c05f39f864e35939c68b))
* update checksum files ([e1cc4be](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e1cc4be0cb154d217152ee45c694f752eb1a4dbc))
* update soft filter to match hard filter ([26dfb06](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/26dfb068c9f036565f033d388f4bd6280d8aaf4d))
* update the pipeline to latest hydra releases ([b668a56](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b668a56aac9b678840a86b9a6bab7e23358e6c26))
* update version for cnv_sv module ([6b52c27](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6b52c272001e15359af6477eb7a1139470870e4f))
* updated tags ([43daaae](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/43daaae1984c936f3be56e49f21dd55aadec9de8))
* use all loh genes ([016978c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/016978c1bdfe5ce056443569836ed1115509cb8e))
* use best segmentation (purecn internal) ([6096cdf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6096cdf640903ac2b721b61418a5846102038bef))
* y-axis zoom control ([450728e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/450728e260a1abf0f034e801a280c80b29a4d190))


### Bug Fixes

* activate all resultfiles ([7b5c082](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7b5c0829ce410da3eaa05d1847289ca10274f299))
* add classes to axes and labels ([ff63504](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ff6350468cf61926a84ea94ea984a58a76a4e6f3))
* add missing function ([7c43eb9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7c43eb915f44284188048df685d56a453800f57e))
* add missing function ([bd0d30d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/bd0d30de6c92d7bc44e719f4153b122f79af5e2b))
* added bgzip ([cd6bdf0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cd6bdf0f58288843b5be5e2aed3751ff367c5a67))
* added cnvkit container for cnvkit_export_seg ([15db054](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/15db054cf70d24bf857e11c826e9d5e4c798ff06))
* added missing rna results output ([d8acf2b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d8acf2be9cdc7bd836cfdb305eb8f6ac9638e92b))
* added tc-file function to input ([9b92fc6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9b92fc674c5fc1d3d9c48e9eb3770e83bd16ef9d))
* adjust margins ([eddf35c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/eddf35c1d512cc6b791718b2153902df5d9236e7))
* bug where 0 became NA in table ([4d68c93](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4d68c934c3c3bd30517a61700ee6a8dae70339fe))
* bugfix ([0193695](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0193695cc461c875410d659c72cb6fd7505f084a))
* bugfix in calculations of high averages ([9e8d158](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9e8d158f05e890a88886622e218d975bdf466e12))
* bugfixes ([59323f9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/59323f9486dfbc27715bd0d2562f6495a99abb1b))
* changed SJ input file to exon skipping ([69e7c4a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/69e7c4ab1ca08fe4cc45ec8b9deff075341f7791))
* **config:** added db field info to purecn ([3be01c4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3be01c44fe3cac35469a4ce1d8b5a0531052c98e))
* corrected germline filter ([6b618c8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6b618c8bf0bd986b79782e981d80b1a6aa848756))
* corrected input file names ([996d913](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/996d913539997b7fec1837c64b0dc875d2f88810))
* don't show table if no vcfs are given ([208ab48](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/208ab48e7675e251c21c779da2dda0382f4fe62f))
* duplicate input files ([2de54b3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2de54b36e4d68837a192c21d6759e4952409063e))
* EGFR exon skipping ([f2e75be](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f2e75bea23f1a8a81d00f46571ad4ebeb2202832))
* empty table not initialised correctly ([f2ea171](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f2ea171084d3f7bc5b28e6496ca131e423e215de))
* file names ([7d4490a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7d4490a8acd7357a660dacae4a6a3c65c05d243f))
* final output filename ([14016fa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/14016fa8fb54790604e15952a15cb36a60a6e81b))
* final output filename ([a6376a8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a6376a805ce6658908f00b615c459a04d0e04c6c))
* gatk_cnv to gatk ([4ca46ac](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4ca46acb5209134acd0559f069402a7d33807a28))
* increased memory for all multi-threaded jobs ([363f0c6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/363f0c6c7d68e2e47ce7afe21126c1ca99496206))
* int + str ([0643343](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0643343b8ffdb79d222b112d34db11036cc77149))
* issue with extra variants being added ([53c65b1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/53c65b1de77b860bce3f7ed177500b43dae4defb))
* log and benchmark ([3c65c88](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3c65c88281e24b94763c7908a02ee0025436772d))
* more ticks ([15e9a9f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/15e9a9f6864bcd88b92373be59405b997cc69297))
* new cnv_sv tag after bugfix ([8f66f68](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8f66f684d495831906b4c7ef8ed013c33bbe48a9))
* new fusion version fixing arriba ([b85d54c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b85d54c920a4593ec6e257e75ac8a5a23d763f84))
* new hydra version ([b904915](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b9049151ffb3b01892ab3cb88334dae2c7a39e3b))
* only include filtered CVNs in JSON ([e83a725](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e83a7253991eadb53bd6473052715ae6ad4c5e4f))
* output files ([d752e69](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d752e69d993f869a8c9f5ea58f254e7eaf829a53))
* pass filtered and unfiltered VCFs to json merge ([fb464c0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fb464c027ad2158f7d5b3161999dbdd51622ab08))
* populate dataset radio buttons dynamically ([eac266b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/eac266be4716c84e1a70b27b1dccccc347d634b6))
* purecn config and updated tag ([c0053cf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c0053cf6c2e28d6f864ccf6e36a08046a9f08e0e))
* purecn rule ([b468cc6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b468cc6c42ac9d4d3912c8f9b19b1693070ba5d3))
* put guides on axis ticks ([fcbeb0b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fcbeb0b407497852f3fe17d54c1114029257bf0a))
* readded modified copy results file rule ([b07020c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b07020ccbdc51471f52232ab708d93cef567d663))
* readded modified copy results file tule ([37804ac](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/37804ace4911dd5f550cdc1e4a7fcbd542cc78ef))
* reference files ([9c7ef2e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9c7ef2e8e09a9166f194e1491a1614e92d0e91e8))
* removed faulty filter ([e9d76e5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e9d76e5c34d865e37a6e318921a0d65081288861))
* removed unused functions ([3b9985f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3b9985ff317480fff968f9f698fc191ff15db9bb))
* rename gatk_cnv to gatk ([194bd2c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/194bd2c10119755ea0b0aca2f6daf2ef5c73fea7))
* revert to old purecn settings ([72acdc9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/72acdc90899c82d6ee0bbd1e80a9cd5d4a2e79ea))
* show what dataset is selected ([c722c6d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c722c6d13cc15a542652491e9f3a82cdc2e03339))
* trigger purecn by adding tc file to input ([3fe8971](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3fe8971bdceed3446572f96e4c7426e2e0891aaa))
* unique copy rule name needed ([1cd3349](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1cd334914f1c2c82310c04030edeecece985d232))
* update annotations version ([5cb38df](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5cb38dfead27214017f026bfe3b0c62f1213cd7f))
* update cnv_sv tag ([fd78dad](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fd78dad781ab415198f5254956bf4cdc658bb475))
* updated cnv_sv tag ([fece910](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fece910d5a07e5efde725347f7b72d7a1b70487c))
* updated cnv_sv tag ([e97edc0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e97edc0fb57af2d6b0fdd55a05412858fe6fc49e))
* updated output directory ([8a52f7e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8a52f7ed2c5ff23775f194c5e753b021a3f7f717))
* updated tag ([4821498](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/48214983c27ba2934a301e68650213afdf506bb7))
* updated tag ([6efe190](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6efe1905d49effd3c02259be6131a232b7325b7e))
* updated tag ([2cd74d4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2cd74d4a6f889d18b906d9bbb0d1a7447c34573a))
* updated tag ([c929935](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c92993577a2196b03af941a53ed20440967f3f3a))
* updated tags ([85e257e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/85e257e354abd3317c9ceb23a79c02e68e7b0c53))
* updated tags ([c0a4afe](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c0a4afe127dd595a8c26cd3ce964bb176c75d09c))
* updated tags ([51191dc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/51191dcb3270057b16e3b7219b4bc02c4e0ff3f6))
* use annotated vcf as input ([44fbab3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/44fbab305527fbb9af0b684d3bc6d99c535e3b1c))
* use blacklisted germline vcf file ([86d4e7c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/86d4e7c0b54e1c294f5addce97217d6e06f057c4))
* use indexed vcf-files to avoid warnings ([b5cb43b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b5cb43b94559e1e0ed428cd1a0ddb9957fc15a67))
* use new tags with updated file names and annotations ([52f9a8b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/52f9a8b8c4c9411d7345aa48fa6672ca04085006))
* variable name ([8700015](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8700015057a9a59034658148eb5d2db457763341))
* writing to closed file ([f180ff4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f180ff49ef5a609aa3af4cad13d3e823d5c9bc28))

## [0.3.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.2.3-validation...v0.3.0) (2022-11-22)


### Features

* add cnvkit html report to final output ([4384b81](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4384b815927e99f3fec2409ea86e7c7ac8f4ed92))
* add GATK data to visualisation ([4407006](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4407006f3a5c35a53f5e04fd178b7f3aae2b0b68))
* add gene annotations ([3f5986e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3f5986eafb9f2083de280da93614d6a77a057f79))
* add gene label background ([2352864](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/23528649acbc4449548d4f3f1083b2fb6c336d9a))
* add report template to schema ([5b46741](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5b467417a5abee37cd27d96b85a35e63ed07a8b6))
* add resource to generate rules ([9593690](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9593690b1a2bc592be93b68dac870b16346facad))
* add table to CNV report ([87cabe3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/87cabe3e5b9f40b83b2704a54a67556d16b88051))
* add three fusion partners to gene_fuse ([fb1f4af](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fb1f4afdf54b0989a25570f4a386800cf3560c08))
* add VAF to chromosome view ([6a96631](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6a966318d6e516b10b27d0c43c88709911aaa795))
* add VAF to JSON ([c0aa776](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c0aa776125da43e9f26e6232599c4e9fae7db290))
* added contamination check to pipeline ([aedb308](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aedb30895f5d86d1949dba2e813efbe22739ca92))
* added fp filtering of one fusion gene ([9322c6d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9322c6d426b32a53effd39227842da0a04e85661))
* added gatk_cnv scarHRD score ([f7df157](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f7df1574c27ad0c50172af3c21ecdf465929f002))
* added merged mutect2 bam files and bai. Also put cnv output in one folder per sample. ([1fab0b7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1fab0b7ad14c55586a94b911a74b74fb6c65a2f1))
* added scarhrd calculations ([f729310](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f7293101017cd42014576dc8caf1f41ae4472080))
* animate dataset transition ([857feec](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/857feec90628ab12a3d8f46f89900cdf2251e91a))
* big data validation for develop and master ([8829aa2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8829aa22ddc382291c52a66fe65cb5bc2e38a4d0))
* bump up cnv_sv version to v0.2.0 ([6c58d5c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6c58d5cd05be0533eb316331220088ecde6edd18))
* **config:** better ranges and coloring ([c4fa71e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c4fa71e293d9ee4c5faaa657f37e7c6b51af95c8))
* **config:** decreased limit for 1p19q ([3ad9cde](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3ad9cde65b1644ea37110e2085dc3cfcaef3bbdc))
* **config:** even better ranges and coloring ([f2abf75](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f2abf758e6817ae4c0585d628b6c7daeab296837))
* d3.js visualisation of all chromosomes ([7fef05f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7fef05f9dff8cbe20f4b52fdeb5c4949130ffc9c))
* enable exclusion of chromosomes ([e2052b8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e2052b8666d574ad46e34316129ba260e7d91f8e))
* **gene_fuse:** gene_fuse fusion filtering is now configurable ([6f7171a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6f7171a9203898d5f3583fea771210fa5991d4b0))
* genome view before chromosome view ([e892aef](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e892aef8008e4a048120a301d8ae7cb99e4bbf3d))
* indicate selected chromosome ([8d74eed](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8d74eed6b900db414c0a4f6882ef25871846c639))
* limit log ratio y-axis to [-2, 2] ([c814fce](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c814fce0c52d2932165bcc0f5560c00743fd3ff1))
* make gene annotations optional ([5ff402c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5ff402cf34757cb77a11b74d7be045ef1463ac95))
* **manta:** run targeted manta ([d38d22e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d38d22e2d63738711e80b9ba4f2ec4d82b4e2635))
* **multiqc:** added contamination as custom data to qc-stats table ([5c5f407](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5c5f40794842655668e1db86e6dfc9664eb5e46a))
* **multiqc:** added contamination to multiqc ([4f7a5cc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4f7a5cc03461e88552ffd61983a7bcc31d68e952))
* navigate to chromosomal region from table ([6e5bb5f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6e5bb5fb54686f300a0c91f078373ea5b43270a5))
* new cnv_sv realease ([d7d1e83](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d7d1e837a6795f2954418930e3041aac90b15fe3))
* new configs for running the pipeline with hg38 ([d4c6fce](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d4c6fce4f886c52c7a5086b9cd139cabcb74d8ab))
* plot all VAF data in genome view ([86e92ee](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/86e92eeb683db8235df4676cdaba6b7fe8ed949f))
* plot selected chromosome ([7f6cff1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7f6cff12903f6420bbab4bdabef1bc44adecbf3e))
* plot single chromosomes ([9d59b84](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9d59b845f97deb5f8a52336c18de617ef0e52887))
* rename files to better reflect functionality ([2200805](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/220080534542ce40a7ea7914f9fee2d84cf1c132))
* rm inhouse HRD output ([58c970b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/58c970b626b1ad6ebfe339f360197957750cc0d8))
* rm scarHRD based on GATK_cnv ([ed8d8ea](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ed8d8ea3b0da8fc715d6c9387351b6b31e0023d5))
* scarHRD + CNVkit only on cnv backbone ([f1fa85c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f1fa85cba4294cabda1bf3d47dcb1c5de2c343be))
* **schema:** set which version of default container that must be used ([633b98f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/633b98fe55723f2407494494e83bc397e5eaa54e))
* **script:** update and correct fusion filtering ([e73dbee](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e73dbeeaaeb15a8574cf778ac8b357b3398770b9))
* set version for misc ([2468043](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2468043e16990eae7e1801e1859658ef0d6e5a99))
* update list of output files ([f47375c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f47375c57d07372f1b0025137d7819b95a01e855))
* update test_input_ALL.tsv ([2468043](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2468043e16990eae7e1801e1859658ef0d6e5a99))
* update test_input_VAL2022.tsv ([2468043](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2468043e16990eae7e1801e1859658ef0d6e5a99))
* updated tag to cnv_sv release ([46f037e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/46f037e0db223a44b133bd6d119a0caabeb2d077))
* updated tools requirements. ([2f0f629](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2f0f62912f84c16f32c81983038ead52fe9e6b71))
* use filtered CNV data ([94a004d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/94a004d1f7cb85bd8feaff033a196d72df75cd8e))
* use svdb vcfs for results table ([fa121d6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fa121d68f0eb3b234255fb621fe404263ae442c5))
* write cnvkit data as json ([009319e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/009319e431a4c9e066d9930f8d0b46a44042aed6))
* zoom in chromosome view ([a6e554c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a6e554ca353078cd505f3dd23b2c019aa401b1a6))


### Bug Fixes

* add conda and container to copy rules ([21da2e1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/21da2e18dc12c637e0e20ee9c004a9253b4752a7))
* add key functions ([66aa201](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/66aa201a4a04cee65c73cd92be326977f11865ef))
* add missing (new) output files ([9593690](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9593690b1a2bc592be93b68dac870b16346facad))
* added folder to activate generic samtools rule ([1e1a59b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1e1a59bb8c6d0282745f404063fb42fb6628b357))
* apply clip path to gene annotations ([b28ec6e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b28ec6e16f68596e5c3e5b386dceb5748e8071f7))
* avoid duplicated genes in JSON ([6d6d433](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6d6d433da32577e24c4df1e244f753a14463f1f3))
* better hover-effect ([5202a82](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5202a82a5d7e8bd018b944e44c2dc5cd3dd56704))
* better transitions ([5e281d8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5e281d89ad077844a7f0fd6bb1ed62e54bec80ac))
* better y scaling ([3d052ec](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3d052ec581c09a14a22b771e7949137d56fa7049))
* bgzip and tabix input vcf files ([cce76f6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cce76f6ef6b9d9913224c67e11808253014e4d6e))
* bug in new snakemake release ([3a40227](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3a40227c48a7dadb54573c26d7de182ba85c2d3f))
* bug where animations were interrupted ([499aaaa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/499aaaa64508a2a3ce8712e39d8fcec3ee0d219f))
* bug where gene annotation got negative width ([af3daa0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/af3daa051bb3da8ea8c7c269f495bbe8e041fc9a))
* **config:** bed files that include all hotspots ([7619819](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7619819bd414d73712773ad5dba84cd7203f0bd5))
* **config:** interval files that include all hotspots ([7725130](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/77251302d701a298a69e0a091b3eb635de366707))
* **config:** rm % in table and corrected column name ([d0598ea](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d0598ea0bfba181bdbd4bc05b8fa65a2c207a06c))
* **config:** update common container to work with cnv visualization ([302f091](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/302f091580452c0f72391df429ae967066096cab))
* **config:** update config to match new cnv_sv release for manta config ([96b2f3b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/96b2f3baf81e6de2643d33ab2518f4d38b8d70f9))
* **config:** updated cnv caller name according to new module release ([6e92714](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6e92714080dcb80b03bf45010f0fb4f80d91b579))
* correct output file name ([5bb11af](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5bb11afa78b4851b4541616a42efc14b7566fadd))
* do not overwrite all output files! ([6413ab9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6413ab9697ce6620ca806806aa389e17cb0a5630))
* flexible sizing of svg elements ([fd74250](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fd7425098c0d5f561ef4093d6dd662b7c1ffab93))
* **jenkins:** keep files on failure ([1be36b9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1be36b9f78c323d0397bbce2bec418b0e0958660))
* less cramped x-axis ([e409bfd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e409bfd8b3b0b3bb4dabd8d0802f60fa4181a405))
* match new hrd values ([4bcaaad](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4bcaaad425702bc982dfeeb7d092ff25c1ce153e))
* match new region filtered manta results ([c224e19](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c224e19f2911abbdc989b0840b989c67d31e85f5))
* multiqc config ([a236093](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a23609388040a1d1fb04bfc7ea4e1e7d9df80a5c))
* new qc tag ([9b28d84](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9b28d842eb0f7599885331a78619ba151c6cf4db))
* new tag for qc ([32ac232](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/32ac232304fd9cee57675d2fb2d64bf4673c789c))
* nicer fonts ([999e78a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/999e78ada0f5ec88ce6bb6706bfb898299cd2fea))
* output file names and updated qc tag ([45fc55a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/45fc55aad378e69c80a8667c20786e9902e1ef1d))
* prevent zooming outside chromosome boundaries ([2063f7f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2063f7f43c9e205d07a8ac7c9493eb68931167b1))
* problem with empty CNV table ([1a0af92](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1a0af929cd5d21c5c45893588bff64685546e34f))
* remember chromosome selection ([511c99f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/511c99f6e96c5022f6b38f57f919e5f08f159e81))
* remember zoom range ([4f18777](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4f1877781847d0b1fc2abc12bb0c37a8e0f2e6cd))
* **req:** added cyvcf2 to requirement to fix cnv report ([cc4a6aa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cc4a6aa15fb4133cfdf04863cfb832727ebeed36))
* rm mgc tag on files ([a1f35ab](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a1f35aba2aebbbd5d7b21485ed8387f0e8bdef79))
* sample name ([478ced7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/478ced7a28da4c913f71b1eaca4ae677ec4f7b4e))
* **schema:** update required field in singularity ([eb48289](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/eb48289bba6c46e8e53e31b80af5ebcc5b8f76ac))
* **schema:** update to match new config ([b60083f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b60083fe05df9c93b31a124e6ac1187148498675))
* **script:** rm extra chr in chromosomes in exon skipping ([c33a0aa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c33a0aa37b524f14d67eed73783ab9e4814359d9))
* set data for chromosome view ([dfb88ed](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/dfb88ed6bf72798ad513284db78a47395f410a5b))
* table styling ([962158a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/962158adf21cc16187eef1977c3a16ee9f193fb5))
* test fixes ([34fc26d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/34fc26d09da1687dfc97c6c101e3a7da54763147))
* test run errors ([3fcc9f4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3fcc9f40408d38e1e3b8db81f3b7445cf6bc4fd6))
* test run fixes ([a64502d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a64502dd566eb995b59b0701fb3fa7e4fb2cee4f))
* update checksum ([7f2708c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7f2708c2f01026341f03d7803aba3f0c049aeb86))
* update config/output_list.json ([2468043](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2468043e16990eae7e1801e1859658ef0d6e5a99))
* update tag ([7e5b548](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7e5b5484c8159bc0d22f5fd4c584f2676a74f048))
* update tag ([152fcf3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/152fcf362775a1a58862b7a14ac549bfaced8c2a))
* x label alignment ([5a6c797](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5a6c797b5d0d190d3da787870fab5556ccace66c))
