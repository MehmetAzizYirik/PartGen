[![DOI](https://zenodo.org/badge/207172204.svg)](https://zenodo.org/badge/latestdoi/207172204)

# PartGen

Copyright 2019 Mehmet Aziz Yirik

## Introduction

PartGen is an integer partitioning based structure generator.

## Method

---

## Download Source Code

It is assumed that users have git on their system and have initialised their local directory. For more information [set-up-git](https://help.github.com/articles/set-up-git/ )

To download PartGen source code:

```
$ git clone https://github.com/MehmetAzizYirik/PartGen.git
```
## Compiling

To compile EQGen, Apache Maven and Java 1.8 (or later) are required.
```
PartGen/$ mvn package
```
This command will create jar file named specifically as "jar-with-dependencies" under target folder.

## Usage

PartGen.jar can be run from command line with the specified arguments. An example command is given below.

```
java -jar PartGen.jar -i C6H6 -v -d C:\Users\UserName\Desktop\partgen
```

The definitions of the arguments are given below:

```
usage: java -jar PartGen.jar -i <arg> [-v] -d <arg>

Generates structures for a given molecular formula. 

 -i,--molecularinfo <arg>   String of atoms with their implicit hydrogen
                            information (required)
 -v,--verbose               Print messages about the duration time of the
                            generator
 -d,--filedir <arg>         Creates and store the output sdf file in the
                            directory (required)

Please report issues at https://github.com/MehmetAzizYirik/PartGen
```

## Running the Tests

For the Generator class, a test class called Test-Generator is built. This test class includes the tests of the main functions. The outputs of the the functions are tested based on the size ( or the length) of the expected output files. 

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/MehmetAzizYirik/PartGen/blob/master/LICENSE) file for details

## Authors

 - Mehmet Aziz Yirik - [MehmetAzizYirik](https://github.com/MehmetAzizYirik)
 
## Acknowledgements
![YourKit](https://camo.githubusercontent.com/97fa03cac759a772255b93c64ab1c9f76a103681/68747470733a2f2f7777772e796f75726b69742e636f6d2f696d616765732f796b6c6f676f2e706e67)

The developer uses YourKit to profile and optimise code.

YourKit supports open source projects with its full-featured Java Profiler. YourKit, LLC is the creator of YourKit Java Profiler and YourKit .NET Profiler, innovative and intelligent tools for profiling Java and .NET applications.

![cdk](https://github.com/MehmetAzizYirik/HMD/blob/master/cdk.png)

This project relies on the Chemistry Development Project (CDK), hosted under [CDK GitHub](http://cdk.github.io/). Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is available though our [Github site](http://cdk.github.io/cdk/).
