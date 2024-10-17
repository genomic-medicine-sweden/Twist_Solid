# Changelog

### [0.15.1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.15.0...v0.15.1) (2024-10-17)


### Bug Fixes

* correction of output text for 1p19q ([e98873b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e98873b37d2e9ea99dfdafb0f6937cb6762adf05))
* update cnv_sv module ([9c1f42d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9c1f42d5a0c3d62c655f5ca761ef6d19a2014cbe))


### Documentation

* update documentaion ([215aa2a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/215aa2af21dae925e38f1ab495b3e279847361f3))

## [0.15.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.14.0...v0.15.0) (2024-10-07)


### Features

* adapt multiqc dna to new mmultiqc version ([a46dfb7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a46dfb796eada3bfb157212d00803eb9f0f26ad9))
* add a rule which adds the FP_FLAG to vcf header ([246f76a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/246f76ad0bab8cdcef2e4563e221759c77314f53))
* add fp-flag column to html report ([de38e72](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/de38e7292ef478b4a94082d35a78ed773746931c))
* add support for changing html report main table ([954b90a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/954b90a1bcecc9fa06de644b3dfe51d746b4c2a1))
* added cnv caller jumble ([13918e2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/13918e2af65de5223f1b72eb98bc7856fa07b1d9))
* added jumble reference to pipeline ([4115c0c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4115c0c0743579c09a9e54b8c633ac9fd851f443))
* added jumble to html report ([a40021d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a40021d030618885ddf565e8e0591b1766e2fe23))
* added jumble to html report ([28a5403](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/28a540383e72d985c4392e411e28d20305c4edfd))
* added matck_cutoff to config ([6405f2c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6405f2c25f26f8f6e51cf3d6ac3653faa2de3ee2))
* added new jumble PoN to reference files and configs ([a986e9c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a986e9c40eb4a1b42dbfd922a2d302e19f4dc7cb))
* AF rna fusion filter possibility ([d46e228](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d46e228b849a7f3f1dfd0546e24cab320c187ffd))
* change calling model ([36dd6aa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/36dd6aa8ad1600b95ac8d4813c25ac31bea89d4a))
* change max size for FP to 15M baser ([0c10799](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0c10799916a02ff8ef8f98d8cdacbaf3780d8090))
* flag FP cnvkit dels and amps based on gatkcnv data ([6ddbd02](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6ddbd02734c07d43276a6e18ff152e8d03a2ea52))
* make max size of cnv configurable ([5fa5908](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5fa5908a324db22eee868b68c75914626f2a7439))
* make sample_mixup_configurable ([a0fdd26](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a0fdd26d0539f16de65d4fb667d5653f641e3ade))
* make sample-mixup configurable ([c6fbf03](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c6fbf03210bb01faae333f457b454bf2444c779a))
* rm CDH1 from reporting of cnv deletions ([32cd1f5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/32cd1f5070ccc2c837952d5f1c0579b4ee7c2b09))
* rm noisy MUC6 from cnv VAF plot ([a15663c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a15663cb60bc462d2d4416a815dc54499f587592))
* stricter BAF and cn criterias for chr aberations ([9600785](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/96007858d634e32debb7c65ac2332ad3e6690d49))
* update qc module to v0.5.0 to get correct samtools stats without bed file ([71a3c0c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/71a3c0c8bafc50188bfbf4829fa74dd6a27267dd))
* updated ENC hotspot files ([3354710](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3354710008e7a9c31b081b70a539e6fe3ca17299))
* updated multiqc for version info ([71e8c8b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/71e8c8b977b0a69a1cd24d6de5e79bc267d06a73))


### Bug Fixes

* add '-' to FP-flag where missing ([6e478f2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6e478f266c8163889fe03dea78568622f049357d))
* add FP_FLAG in header of all cnv vcf files ([adeab23](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/adeab23e265dadd6505619904ccfbbcd66aa8697))
* add missing python module ([8566c04](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8566c04bac21cd61b0ba37849a4c305712c590ad))
* bu ([ff962a4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ff962a4c5ebe9415915368920ba5416d2f705e05))
* bugfix ([efa90d3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/efa90d39ae4015a9516f348e5f63b168faf0b251))
* bugfix of version crash when pipeline runs over several days ([2407c7e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2407c7e71ab8cf6f85e13cd6a81996ae3df3d6fa))
* bugfixes and improvements ([aa84845](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aa84845519c0464ac6cd037f71a3022993747bd0))
* correct git version ([7bd5582](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7bd5582245a5dfc989abba88d9677c315dcde88b))
* correct input to html report ([309667a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/309667a34a9c5405f1074bb230ffcf941d9e13d1))
* correct output variable ([8b5c13f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8b5c13f1ecb4c3620cbb3a0e53f03ee5d53c0804))
* correct reference version ([8adec28](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8adec281403e9d6c2d7378f789da9d9d7ef85d61))
* correct style of table header ([0edd510](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0edd510b8e8efa22b8d035422d10a093bb25c404))
* correct style of table header ([62cd31c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/62cd31ce4b7e83983d7a0e132deed35b6389fecf))
* correct test output file name ([a811096](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a8110960703a446308f7ef0d384b70744ead37f5))
* corrected cytoband config parameter ([94678e6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/94678e60f7ba217dbca8af61f3467ad9f3815203))
* corrected git version ([e6c8da1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e6c8da111f2f754e1d2a88d5c4c0f4e859f4cfb3))
* count files in correct folder ([0c700a9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0c700a99cef1d3c1b5e33931236321be767771ff))
* do not use ? in flag ([994de7b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/994de7b5433758a82f4d4eead7c245b1c59f23ab))
* file list bugfix ([d8cffa4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d8cffa4abe6cb8bc3f08a626cb381068f53bde3c))
* final fixes ([a75a249](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a75a2498a6432ed09ee634b0a0b86e2dbe9c1765))
* flag all variants with - ([220324e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/220324ed775b192218f3f75e64e1eb8ee2e9ebc3))
* fp-tag PR mismatch ([0f106d1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0f106d1b884dd0d03237a1c37b3f0df4351d6092))
* FusionCatcher liftover fail crashes report ([99a9e34](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/99a9e3458f177d82bf3e2df7bb914f63bd16e9b7))
* handle missing ID-SNPs ([a247ceb](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a247ceb1b71bb721f6a6d98a0a63b3d90ac31169))
* ignore version files ([d3d0e42](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d3d0e425a15dd2a36fcae9e3c546b3db0e1471de))
* index vcf in for pysam ([8735f9d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8735f9d7ccb024ea85ca3f799bc4e12d27066e78))
* jumble pon override ([c261b88](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c261b886a79a209099c5f1d617a89f125566354f))
* multiqc report adjustments ([e05c06c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e05c06c9b912fe278e3d1125713dd75147018649))
* pysam bug ([1927d0a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1927d0a3f4ee528d8dac966b806c6996ddf0b574))
* pysam error ([c4331f0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c4331f08f9669e1a1d532b6337702659fe4a6222))
* ref files and configurable AF filter ([dbf208f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/dbf208f1fe15de818321444bd8a7ef7d8f3cb0cc))
* remove duplicate input ([206715a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/206715a6d93cec7f1a2034cf967a32e32045f3ab))
* rm insert size from rna multiqd ([bdeebd9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/bdeebd91d825bb9321eecb10a1ca2026bf01acc7))
* small changes ([984db35](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/984db35d29dc460cf85c258401c81f239000d86f))
* solving the FP-flag bug ([2048307](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2048307866b2e89714b792c81d4134a7204352e0))
* update annotation module ([04b5232](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/04b52326f2c78233b494f853ad20d7707047970e))
* update cnv_html_report rule ([8abad6a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8abad6a0db60034e5169139e478cab78485ffe04))
* update cnvkit rule in Snakemake to new cnv_sv module ([74a86ae](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/74a86ae080456ab46dc5ae2487d21543e6aa6c00))
* Update config.yaml ([befb670](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/befb67058d888b0be62c3817f4fd0c8024efe7cd))
* update hydra min version ([b924f43](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b924f43845598d80b36365150a19aec9f5ea2757))
* Update output_reference_files.yaml ([b2e8c6e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b2e8c6e11a03f833779c033252356563037934f3))
* Update resources.yaml ([62e766d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/62e766d993531f8c3bf72c8e252e4ef080538d23))
* update rna multiqc report for new version of multiqc ([991f187](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/991f187e13dd9ef4485656add990998982cbfe4a))
* Update workflow/scripts/sample_mixup_check.py ([c121a4b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c121a4beccaf29ef9fa5b11bffa9b404c9d57dca))
* Update workflow/scripts/sample_mixup_check.py ([3ab78f3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3ab78f3d1e6389a294f0ec7ea23319752aa8d890))
* use samtools stats instead ([535d53c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/535d53cf34dba432fdaaca2161e03fe7640ca742))
* wrong fp thresholds ([c223d50](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c223d5014c9c94b6f797bdbaafcf786a4f088b6c))


### Documentation

* corrected purecn documentation ([8a4ed41](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8a4ed41a9a34c963ddcf5404d9a86ad3fc0eaf70))
* sample_mixup ([9527682](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/952768223e4624b4340d6c90767385320ed3ac44))

## [0.14.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.13.0...v0.14.0) (2024-06-11)


### Features

* add hotspot file for ENC ([2de88ed](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2de88ed957e93287d806fa31850ac55baf4ce5bf))
* added blacklist filtering of small cnv deletions ([f6e77bc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f6e77bce995bbdafb2ab022e008c2d98e1a2ab6c))
* added ENC hotspot files for hg38 ([0e06818](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0e068186e8bad091ddbacd37f24511cf81f94d4d))
* added separate filtering of fuseq wes for umi ([8b83894](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8b83894f01aca19ef1e679c4353b336164caa2d3))
* new PoN on figshare ([41f3f0b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/41f3f0badc34b9200c43c640543776df9037a217))
* new PoN on figshare ([f702d7c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f702d7cf3e580455456599fb810adac42427509f))
* new PoN on figshare hg38 ([641e37d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/641e37dc71de4b83d5a1d3d3736dd41e74f99edf))
* reinstated genefuse for umi only ([a5af752](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a5af752479485b6bafae0d3dbf91fd211fbd120f))
* update common container to new hydra ([6f64d68](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6f64d684eebc222feb6d8de0d2ecf9188a01b9f8))
* update hydra version ([651c38f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/651c38f2bc8ab83fc426220dd534e3e1a2799084))
* use latestest version ([365d3d2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/365d3d29f292bf3c7e402c68f3fa17fc979acce9))
* use new PoN ([e5bbd55](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e5bbd55f04d40cdaddd9474b234b0e0d2eddf1b1))
* use new PoN hg38 ([1151f3b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1151f3bbf7b8b5043b2f3e72b2f68e4528a7edf3))


### Bug Fixes

* add missing singularity ([24eeac6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/24eeac695380da19ec21cf78558a6b58bbaed6d0))
* added blacklist files to reference config ([3450166](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3450166a98467e10c4bd26e1cd12f2dd563f44f9))
* config schema ([e727bd5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e727bd5579a5d8c325438d8a02a15aecb718543b))
* config schema ([9470d56](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9470d56e06f324c653ad5b096a143477a3182bdc))
* correct reference files ([4ff7a76](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4ff7a7605ba322dfd20850b3386dbb298a61899a))
* correct URL ([ee7c8f6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ee7c8f6d2d3d2fc575bcd83df07d8f7d84a06a77))
* corrected config ([20fb88e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/20fb88e55d34d9f7e11bd6202df087fdf17420d8))
* corrected novaseq hrd PoN download link ([6ea893a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6ea893a407f26fe0ddbcddffa4936c055489f3af))
* git version ([4efd11f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4efd11f3c1ed331244674ce034079ecd2d47ee77))
* git version ([ef0891d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ef0891dda5f7f7c48137186895438f41c0683fb7))
* rm conda ([210dfdc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/210dfdcd40c82d9b41221770bc7992249a5b4d80))
* update ENC hotspots ([fec13d5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/fec13d5326682753b8a5a5465d7b15a553f79b6d))
* update ENC hotspots ([8334f7a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8334f7a2f1e35d91143228cf95e64c4c1dde5c26))
* update ENC hotspots ([66dd017](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/66dd017e9b7348f1d80d116359ea62790a3c8bff))
* update ENC hotspots ([1338508](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/13385087ceb9169e1fdd5cf60e35591aa4490dfc))
* wrong release of design files ([cc246f7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cc246f74d16b38b1ef719a51c7dd1b53face46d6))


### Documentation

* added documentation ([87fbda9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/87fbda9103cd378da03e1c092d700d3816b7d3ec))
* do not activate python env ([8d28cdf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8d28cdf09a512b016551820188a68fc405de58d5))

## [0.13.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.12.0...v0.13.0) (2024-04-26)


### Features

* added calculation of cnv in chrom arms ([667d53f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/667d53fe56eb15047ca326a91957bbfdb30065fd))
* added chromosome arm support to rule ([ccbe735](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ccbe735010ffb680a5d5ae05a1d848453916ee8d))
* added configuration of polyploidy and baseline limits ([8f91b42](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8f91b42a5509fec983c1c077d53c499c035e955a))
* added more config variables ([9a92bf2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9a92bf222bfb160804748166d49b3d54c6b646ec))
* added table to cnv_html_report ([e577cdd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e577cdd192036e240570d33cbc2a0d9cac8f59d3))
* added warnings for baseline and polyploidy ([9a25cdc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9a25cdc0843f188fe1725231740db65d21181d0e))
* decreased chr arm fraction ([b1ce417](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b1ce41711b33a9e076b7e23c73698d36f3244bb0))
* do not merge regions that are not identical ([1d9a027](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1d9a027bc274a475bd054f5f180c5ada3fe8a1f1))
* new amp genes for hg38 ([c6f7b5e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c6f7b5e2ef6aafd0d7dafcd50e1d6330cc8affc8))
* two new amplification genes reported ([c8264b2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c8264b2f887da62bc7231b580702b935a49a4a91))
* update ref files for hg38 ([0ad022c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0ad022c8f2513a5e0becbe967aff0c2d11c2a86f))
* update test file for develop ([36bab94](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/36bab94d269fc31a756e2f211de78afe512173a0))
* updated with new reference files ([5e057e7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5e057e78889ce3e6d331704ca123e6fb660c85d1))


### Bug Fixes

* add "_" to regexp so that indels are reported correctly ([7549036](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/75490367490628b403907cf1e5126483a1435fa2))
* add implied msi configuration from module to config ([a9dea4e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a9dea4efea6fe38490bbd3b656624735a8f07341))
* bugfix ([e681fee](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e681fee009e69cd73b8e78c05e5aa9fc6cfddd27))
* correct missing part of vep ref file path ([36b3b0a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/36b3b0a1eb61ec8416cd90601c08ae957a44e22d))
* handle CNVs bridging the chromosome arms ([8a508ad](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8a508ad13188648eeb08f4d38d5870cf246aa8f8))
* remove unused library that break requirements ([86afc7d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/86afc7da3413906a00c538cdad86372a18e28a51))
* rm duplicate entries of cnv in report ([3a2cae4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3a2cae4c773f5f8d4ec5ec28f1bfc8f1fcb6c723))
* spelling error ([46aef96](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/46aef96bf47f0f5d17b5c18b8994a49c1cd905c4))
* update config, merge_cnv_json, to use latest amp gene bed file ([be98e22](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/be98e22671d1f31d684d4ea45294d7b89f0aa55f))
* update hg38 ref files, point to correct release ([e8efb8c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e8efb8c2162cbb71df1c3ab01ed7c26e28e411a9))
* variable baf bugfix ([d460172](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d46017212ef9a7cc077e090bb11d3c56c9dfb4d3))


### Documentation

* added config schemas for new params ([62c832d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/62c832d0979382d82233eb59421ab5ee6dbf7e14))
* added schema ([b1c9bfe](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b1c9bfea96a960aa0e048753c27fa9858e5c625f))

## [0.12.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.11.0...v0.12.0) (2024-03-05)


### Features

* adapted configs to hg38 ([7c601ab](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7c601abd0a4ec6634b53b5711be7c314353ac1e9))
* added new rule for rna dna mixup check ([164baa6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/164baa603d4026b7660187e15f78bf0777752a1f))
* cnv html report additional table only includes additional calls ([e3971eb](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e3971eb054e0865447c2139ec8f125689d50b030))
* data config for HG38 (missing HRD, SVDB, purecn) ([7376033](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7376033c86cf8841bd75de19a8a2cc047e831f30))
* increased read support needed for fuseq_wes ([f669bc2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f669bc26a40ea3ae11dfc4f64b957a3ff8b97559))
* print pipeline and software version and config to results ([d643a39](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d643a393387cb975b6537da85e0364b45bb372de))
* put mixup report directly under result as it is both rna and dna ([c9cdbc5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c9cdbc5e21a8ace20b77dbb664bb47fe3f78d939))
* removed gene fuse from config ([c1f25ba](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c1f25ba5beffcffc810c0838880ae861ca94d6e4))
* removed gene fuse from pipeline ([83260f6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/83260f6a376abb8d11bf7b0be89f2a4f2d433bcf))
* update hydra-genetics version ([c57d3f3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c57d3f3da9ead5f9868243d36507708605506c21))
* update reports module to v0.4.1 ([#398](https://www.github.com/genomic-medicine-sweden/Twist_Solid/issues/398)) ([32e40e9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/32e40e9bf8916b49f658d2e8ef65f533890a503d))


### Bug Fixes

* avoid clash with earlier filter annotation ([9d2ed1b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9d2ed1b103b3c41a4fd54225c8a247ea52ba322a))
* bug in Snakefile ([d8307f4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d8307f4c1b38a54510aad1ea3a4e788a2167e910))
* bugfix ([a154b11](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a154b110c5b78768f492e10c84ac81ccb2969181))
* bugfix ([75d21be](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/75d21be979cb0803def7a54e6f559ecf02891f4a))
* bugfix and improvements ([dc7f4b4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/dc7f4b4202b2762a17cc7b15daa8033054ec64eb))
* changed sample mixup från txt to tsv ([5b91e60](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5b91e6068fca5ec776f745539c8888f3c180a56c))
* correct paths for CNV report extra tables ([#401](https://www.github.com/genomic-medicine-sweden/Twist_Solid/issues/401)) ([37a461e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/37a461ee645b241405009854e8f2cccd12904ca7))
* extra table without tag to cnv html report ([765ce51](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/765ce518879b3a9cb16938d34a47c47771bf0cd2))
* handle missing files better ([a0be997](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a0be99720cce840688fd701c664d0f50fd5504af))
* make umi tmb similar to tmb ([c525c0e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c525c0e78fbafb80a1acc308767f3e5042b52d4a))
* min 1% AF ([ef2bcaf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ef2bcaf3cf842166d1c779f6800712e89e8e06b7))
* missing vep config in integration ([e8e8a57](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e8e8a5749df2b7fef449d5468f21704ca4e65b35))
* new smart_open release crash ([03e0a0e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/03e0a0e7340c3936cedbe1c040f9d9315f04d890))
* pycodestyle and unittest update ([79cf8f8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/79cf8f8267f425ff3aff12643d9476dcdad5ba15))
* required variable with corrected indentation ([115e018](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/115e018027cfa2a2366cdb6b1068b1d94e424ee2))
* resolve conflicts ([6343bd2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6343bd236a96841b11c5b9835107cb48ad40c533))
* rm dangerous qual filter form umi ([c9d4584](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c9d45840a4c4efcffda44e0654c157571b6558cf))
* rm dangerous qual filter form umi ([d20f0cc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d20f0ccb7b3d9499650a75d5d17aad0f9cb0b6b2))
* some hg19 file leftovers fixed ([5433b16](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5433b16de07cc0e9d865f9ed5e98ecfceddabe1a))
* VEP cache move to config.data ([a36b570](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a36b570d8e62d31691ca56a6b6c62dbfabd8ede1))


### Documentation

* added missing rule ([a6712f6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a6712f6a66333920c756399036a73e0eb218e1d4))
* fix documentation and codestyle ([7d793b6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7d793b64040d3b761d2d36ed21147b9cd6b71c35))
* rm gene_fuse from schemas and updated documentation ([a991a0d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a991a0d85b0c593a34ca36f348a62df08cdef136))
* rule plugin version updated ([2d63ac6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2d63ac6c3175cadd603bd9a1a96a0d2e9ddde69f))
* updated rule documentation ([43cfb5c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/43cfb5c678935f816de60302429fdbd3450a84a4))

## [0.11.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.10.0...v0.11.0) (2024-01-26)


### Features

* add amino acid change to coverage_and_mutations ([a66f501](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a66f501177ac9393c452ba8b49c8202b69b1dca2))
* added dedup coverage ([01046e6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/01046e6841ce62457041d12538c645a0fc40073f))
* added umi vcf filters ([dd1cbe2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/dd1cbe21ba1dd93218771c62500d7a7123962bc6))
* decrease germline support needed to improve cnv visualization of deletions ([c6cd922](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c6cd922d97e4c1899ceb3dc88cddb55283d8e9eb))
* decrease support needed for germline SNVs ([efb2827](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/efb2827002a48bf9a6b94796b4f5c472e02a9afb))
* restructured rna fusion report for improved interpretation ([373021c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/373021cd2b49dfa43103b7f34bd0d3c669b0c010))
* update Jenkins main test ([5b30a35](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5b30a352a9951edda2afb2c5a3cb812e62fd9bba))
* update reports module to 0.3.1 ([#389](https://www.github.com/genomic-medicine-sweden/Twist_Solid/issues/389)) ([9346b03](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9346b033c722a797a67a1f39ec069857f57948ad))
* update to latest common container ([eef00aa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/eef00aa9729bfde0f2378835b1efa25e412a1312))
* updated hydra-genetics version ([b7eadaa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b7eadaa05f046f957838870c47cb4da658b2ce35))


### Bug Fixes

* added missing double mutations to report ([ed1a526](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ed1a526d31713c12bb72114c7a6f242d8fc32b4f))
* bugfix ([e9f2511](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e9f25111a78fdc576b718daa80206216cc3d4654))
* bugfix ([2dbc810](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2dbc81070d027003beb2364b53c4529df93733e2))
* bugfixes ([7fcd240](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7fcd2403b681e3ae62076e52407903d7a6f029c5))
* changed sorting and bugfixes ([355d99e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/355d99efcb0d823996b67a22ba346d6cddaad15b))
* corrected filepath for rna in develop ([85abf47](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/85abf47bf827592b47ef70cc0bf08017d70c3993))
* corrected input file for hotspot_report.smk ([8ce6512](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8ce6512d526ce2ffcc11c71055654b5dcbb9e002))
* corrected input variables for hotspot_report.py ([076f5fd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/076f5fdb35ce124d36a69e2244f086851e36bc03))
* corrected input variables for hotspot_report.smk ([6758fa4](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6758fa42d0b05d192713ab0269e5f2d617df036d))
* corrected results path for soft-filtered umi vcf ([a8ed44c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a8ed44c90a8dca842d1ba7d3fd837586f6d413f4))
* delete extra file ([2803f67](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2803f67ef1141d6ce224f550111411ea62fb1324))
* housekeeping dict ([7105a88](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7105a8830028f0a57578f01be4203fee22d16a78))
* input file bugfix ([bfa94a6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/bfa94a6fd1bb09a9b17b30b463307d4878409da7))
* pulp version ([f37c08a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f37c08a42bd43286d74e23985bd738e18b9d26cd))
* remove extra header columns ([6bb5d42](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6bb5d42afd1dc4f3f702302b75079e73fce73e98))
* remove extra sample folder in output ([39fefd8](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/39fefd85154a62b31c7c41654d5abbc964f46b15))
* remove newline after purity in cnvkit_call command ([b820ab5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b820ab51ec10b8685cbe36afc482f9146e5174c5))
* remove unused time consuming step in fusioncatcher ([635a40b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/635a40bcf62b51dff051322c40bf41e127cc0e8b))
* removed spaces in jenkins test develop ([05eb8ef](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/05eb8ef564a97f80c7a386c6ef545aadaac50818))
* reverse sorting and add chr to fusioncatcher ([794f8f7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/794f8f7c6924c019d5139cf38b3de9b6dcb287fb))
* update hydra-genetics to v1.10.1 ([be2682d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/be2682d79c8c37a6f43be7940ff4d3166e3e8979))
* update jenkins start scripts ([8cdaf1f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8cdaf1f3e0500e1c4af308bd9a485f4c6d6b39fe))
* update profile with mounted home ([7665c8b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7665c8b4108c24558d030304ca3f219008969fc5))
* update profile with mounted home ([0b1be40](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0b1be406eefd3d36ea3324dffcfa98e941847194))
* update profile with mounted home ([1947e3a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1947e3a8d170250146529847df7dac2e5f6d884e))
* update profile with mounted home ([b3d56e2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b3d56e2fe3a2227789b947889ce4ef818ae1b749))
* updated jenkins output for develop ([887584b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/887584bf8775e15a54f053b62d93b1bf9bd802db))
* updated to bugfixed common container ([1449abd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1449abd6b6e8b74623985bd7b22073a225346cf7))

## [0.10.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.9.0...v0.10.0) (2023-12-01)


### Features

* 2 new TERT positions added to report ([3617acb](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3617acb951c9d471cb3903fb536982c4aca01733))
* add blacklist for fuseq_wes filtering ([6bc03e3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6bc03e308bf08b16dca0d5cbcba72f8e83066c3c))
* make fusion filtering configurable and add more fp fusions ([5ed5c97](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5ed5c97c032e28027f384c095bb867368958fa8c))
* more conda updates ([df2b77c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/df2b77c18fa3aa38c8f270fff3a9e30bf4aee3a8))
* report internal callers for fusioncatcher ([d5fccad](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d5fccad6e6a908274d1d9d72601ff7b399bc1922))
* setup files for reference file validation and fetching using yaml files ([7078c04](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7078c043a3e05310d70a898b017d30f66df10d6f))
* singularity build script update ([8a8b665](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8a8b665e178b4980b7c6ab6dffd6db6fbba92f86))
* update to latest files used by config ([07e351c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/07e351ce02df7a2cf623c19ce305c0497d60ac7f))
* update to newer hotspot file and clean config ([72ef4ff](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/72ef4ff102c3eda03e2f175604994aa0b2abc942))
* use v0.1.0 of fusion module to filter FP Fuseq_WES calls ([3c30fb1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3c30fb10cde2b23a5cfe86a56236c8d123e59621))


### Bug Fixes

* add lambda to functions so that they are only valuated when used ([0b85dda](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0b85dda077d46451074b3cd34d8bd3b8e3e3a5ca))
* add missing extra option in star-fusion ([7476a0d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7476a0ddcbaf72d6aa040e989b291e37a8e30f67))
* add missing type ([35e5f6e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/35e5f6eba600695d52de07bf3f9b7f6c7df611ac))
* change from svdb_query file to svdb_merge file ([1ed29ad](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/1ed29adbd2474204bd71eed92d95215faa89d502))
* configs for reference profile ([6262775](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/626277567977673149cf341ca51ac279df5889ca))
* correct path in config ([dd63afa](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/dd63afae10de654765bea96f29375059f7632a50))
* handle reference creation with unit file containing both T and N ([61b7764](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/61b7764bd502cd28d47d631effb90fa11a5b955d))
* increased star_fusions wall time ([5ff15de](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5ff15de1e6763eda2057650f97ab544f0e56bcdc))
* lamda in correct position ([128aa8f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/128aa8fd206e20a9f9828ea61b634a8a97b56fa0))
* minor config updates ([8296226](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8296226be5383ef1fc887f83ae51c26895caa71f))
* new rna fusion filter file without duplicate entry and header ([280fd4a](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/280fd4a59afb54952ab4b8420d37d1827bba50ad))
* pycodestyle ([36eed32](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/36eed3242ee439363ee7249320a8433a6ef7b9e7))
* support header with fusions filtering ([c8a8bb7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c8a8bb734614c4686ac3251ec71bd7a2eae4d283))
* updarw config files to match current setup ([0e689db](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0e689db20ee5e19ce3746d554a5f2b0f49be2ccc))
* update configs and restructure ([7f397de](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/7f397de7c49d78d606d229d8413d0d6c3dd2dec6))


### Documentation

* update documentation ([0fd09ea](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0fd09ea58561c38466cbdf1cb97243ec7eb0704c))
* update reference page with instructions to download reference files ([3abac1f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3abac1f9fde475163371fa6e6160ca071b2195d2))

## [0.9.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.8.0...v0.9.0) (2023-11-01)


### Features

* 2 new TERT promoter variants ([ee72af1](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ee72af10749fe1f770243e91850291c052f7d9f9))
* add umi support for cnvkit ([979c816](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/979c816a60040871642e44d172ff95306c488d1a))
* add umi support for fuseq wes ([a4d7139](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a4d7139a37dfe5254208d37e1c9b20509b721b70))
* add umi support for fuseq wes ([522ac37](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/522ac37a144a04566967308e7e0e14e691f43dca))
* add umi support for gatkcnv ([a2ce975](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/a2ce975ed69e384af67bc4f3d00007f90b06a3dc))
* add umi support for hrd and msi ([20560a2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/20560a286dc369f22c0afecb69a2f58f499c40ca))
* add umi support for manta ([8a703d9](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8a703d9cd60d538bdf170dd7b6e30519215f5a5f))
* add umi support for qc ([89f673b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/89f673b4ad2aefddfd844996fc80b2c4dd4eee88))
* add umi support for tmb ([9683d6b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9683d6b16eb9f1e62eafe149f8051cef2d41e54f))
* added tmb gene filter. TMB uses hard filtered file for input ([9e242a3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9e242a38e987b4f9107781f898b64bb1270fd527))
* added umi choice to the pipeline ([8e13b38](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8e13b38bbb5ac74204d33f59210c08a0f9c24c4c))
* changes hard filtering to let more variants through ([70606f3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/70606f3e348a2c00b977e66658c14df87fe8a74e))
* hard filter for qci ([9ae0a2f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9ae0a2f0067e76951a2ba13c4b46cbe9820ddc2d))
* min vaf configurable in vardict ([303b5f7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/303b5f7abb99085fd8c75646b541409e84be61c2))
* rm need for ruleorder for copy rules using global wildcard constraints ([ad15f6c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ad15f6c4a7e42c6aea2da4c0ba87763335aca596))
* run msi w and wo umi ([6766b18](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6766b187a9f7cfbcbf949358208216b7e9c8b82d))
* umi in rna ([6ccba4c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6ccba4c0a011c9db6e7fcf62b77a0bd710531427))
* umi vcf filtering based on sample.tsv ([4fdd999](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4fdd999890a55d1f3dba50b0bee63638ff6bbb1a))
* update alignment module tag ([ed33d45](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/ed33d459affe04e9ec69117979d5724f3bdc6f4f))
* update to v1.0.0 relaese of annotation ([3fd0c7b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3fd0c7bd56b70e88952093636e4e7d7265c7f29e))
* Update workflow/Snakefile with new alignment tag ([e74398e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e74398e731c9a9d5493ea8d27c9387a3b78ace4b))
* updated alignment module tag ([62500dd](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/62500dd79cb6f4fad6627b93fa333d69b458e9b4))
* updated alignment module tag ([341a3ba](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/341a3ba3ffae977bce862fbc2f8127955ee236bb))
* updated module versions and adapted to these new versions ([c4aa705](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/c4aa705018c8d533e906d7fee31adfb3d7b5e2e3))


### Bug Fixes

* adapt to new umi alignment module ([60a90b3](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/60a90b38cc243baf289550114c81e173716237b5))
* adopt to breaking change in vep rule ([96b4b26](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/96b4b267fc0829f4b15ffc14228f915f537034b3))
* bugfix in get_vardict_min_af ([cd503b7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cd503b74bcbb91361d15c843b7e8bbaa0f3fa1c3))
* bugfix in get_vardict_min_af ([44d4b59](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/44d4b5958b2a49b4041e80a9f6e113e29b0c8ba0))
* corrected gatk_mutect2 input files ([442d62c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/442d62c93c013ae49eb6f773210eb5f40b9bbb22))
* fix manta output files and rm unneaded umi rules in Snakefile ([5ec7a33](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/5ec7a33699c34fedc03556ceac549e1bc14c24f2))
* fixed units.tsv and adaption of reference pipeline to new annotation module ([9d0c926](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9d0c926319d6d19e8a2b8f4c65c09aa6a00a995b))
* gvcfs now have mosdepth umi coverage ([cf3e9c6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/cf3e9c69abb3d1112bc28975e4bc8132c042638b))
* manta to use original name, filtering wo <=, msi-sensor w correct output name ([0a9960b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0a9960b4506e7920e4e940c2d3e5e1123fc8d909))
* match output result from reference module ([e2415d5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e2415d54bbf204516555be2a2b862c665e3d8431))
* mosdepth should be run on entire bam file ([8643bb2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8643bb20a72cb7ff78589cab2115a9549116dc8f))
* msi ([6578e7d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6578e7d4a0ffd6dfdb6d7ada6d0b7b21026b9a35))
* new common container compatible with filtering module v0.3.0 ([28256f0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/28256f00c5f8d049e235e4f5e62fd9d64b04d172))
* reinstate output vcf file ([2b6607b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2b6607b20bd488bded16266faa68630d439621dc))
* rm bai input for picard insert size metrics ([4f8915f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4f8915f4d233a5c021560be0e0d1a5b7bc05b1ea))
* rm umi on rna ([086f8dc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/086f8dcdc3bef4456d23517860a34d69ce65b002))
* umi fixes ([9ce1a0b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/9ce1a0b8b4f5fa32f30fbf8aeb2462e0cfddbe61))
* update tag and only use umi for snv and indels ([d6b1edc](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d6b1edc7f2d64bdb0dfb7658dc2244a3e885deeb))
* updated alignment module tag ([6207ed2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6207ed22e64a1b29a50f48ef790460ea576c699b))
* updated vcf filters ([def213d](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/def213d7acb25e46f3e6e739e4d51b95aabfa437))
* updated versions of snv_indels and annotation to fix compatilibilty issues ([974dd27](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/974dd27822368d149a155a96e150ac591f075492))
* use same hydra as in common container (v1.8.1) ([d3f401e](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d3f401e20f1750842fbfa222dfa2476388244144))

## [0.8.0](https://www.github.com/genomic-medicine-sweden/Twist_Solid/compare/v0.7.0...v0.8.0) (2023-09-22)


### Features

* tmb two noisy genes filtered ([45397d5](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/45397d508ab84f11d091fd039f1115cde18691cc))
* trim reads to 100bp for improved results in Arriba ([aa06413](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aa06413befbc812bae78e0787a79128350853245))
* update biomarker module version ([d9f0059](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/d9f0059a76bc13fe3c3a25c16810227386e78437))
* update reports module to v0.2.0 ([3f9be85](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/3f9be856ce2c8cf415f01b0606d70c0687a858e0))


### Bug Fixes

* added modified input function to cnv_html_report to fix notemp bug ([b535c42](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b535c42c62914460a32b73797eb521d3c86380df))
* keep bai files ([753b789](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/753b789d3ff2d4aec8889903a7e9b26f8b427d31))
* merge two rule definitions into one ([e294971](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/e2949715dbbb33ae8b4d077a34bddc5b1f7f6196))
* update snakemake version to avoid checkpoint restart job bug ([67f1bf6](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/67f1bf6dd693078d47de9745759e485bcd502239))


### Documentation

* correct rtd links ([46de326](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/46de326494df35c6e6ab004a628b9106b846f958))
* correct rtd links ([f793098](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/f793098c03065b3595b2f264e8f1e44464eb549a))
* correct rtd links ([51c0a55](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/51c0a55580352a2ff8a642e91a69f0e203c79d0e))
* correct rtd links ([8e6a3be](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/8e6a3be0c244d19a99c621a8687f3d2fae245460))
* correct rtd links ([216e0b2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/216e0b29540aa7cc21d559faa0d8ad4036bab14c))
* correct rtd links ([20b021c](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/20b021c67342a72542c66618296a20458f7d810c))
* correct rtd links ([497db09](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/497db093bd8ec7378e834378f68736bea05db5e1))
* update PoN description ([92db1ee](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/92db1ee65d24e22a865112ed56a640412747fffe))
* update readthedocs links ([aac222f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/aac222fccd3e18882f885c8c1baf50362e0729c9))
* update rtd links ([0347a62](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/0347a6206d161727e12058ef3c8cb139889997bb))
* update rtd links ([576112f](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/576112fe1417b9e210ca1c1bf993e5a849c1b372))
* update rtd links ([4f07de7](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/4f07de7454989cf197f5a1426c215012bb377d14))
* update rtd links ([6f99f6b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6f99f6beb226765f33a7bf96dc3198740d79670a))
* update rtd links ([2c353c2](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/2c353c259a013585746d0eadf5031c1cc48b7e4d))
* update rtd links ([6a1c7bf](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/6a1c7bffbe82ef2f64412a860ad183d01440350f))
* update rtd links ([238d74b](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/238d74b02fe62a2ba033345721c1c9a78bebd409))
* update rtd links ([80b0637](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/80b0637e40d893153f14c90f994bda40256a9f17))
* update rtd links ([b2e3967](https://www.github.com/genomic-medicine-sweden/Twist_Solid/commit/b2e3967cf99579c119cf7b294575f22c93d36769))

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
