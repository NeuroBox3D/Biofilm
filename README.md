# Biofilm Plugin
### For support or additional information feel free to contact Bettina Buttaro (bbuttaro@temple.edu) and Gillian Queisser (gillian.queisser@temple.edu)

Bead-Evaluator: [Bead-Evaluator Project](https://github.com/NeuroBox3D/Biofilm/releases/download/v1.0.0/bead-evaluator.vrlp) Biofilm Plugin:   [Download Plugin](https://github.com/NeuroBox3D/Biofilm/releases/download/v1.0.0/biofilm-plugin.jar)

This repository contains the source code of the VRL-Studio plugin developed for the publication [Development of a new bead movement based computational framework shows bacterial amyloid curli reduces bead mobility in biofilms](https://jb.asm.org/content/early/2020/06/23/JB.00253-20/article-info).

### Using the Bead-Evaluator:
- Install [VRL-Studio](https://vrl-studio.mihosoft.eu)
- Open the [bead-evaluator.vrlp](https://github.com/NeuroBox3D/Biofilm/releases/download/v1.0.0/bead-evaluator.vrlp) project in VRL-Studio
- Set paths and properties and invoke the desired computation

### Download the slides from the tutorial [HERE.](https://github.com/NeuroBox3D/Biofilm/blob/master/help/resources/img/Bead%20Evaluator%20Biofilm%20Plugin%20Tutorial.pptx)

### Video Tutorial

[![Screenshot](https://github.com/NeuroBox3D/Biofilm/blob/master/help/resources/img/Biofilm-play-tutorial.png)](https://ensemble.temple.edu/hapi/v1/contents/ebb438f5-be08-417c-b830-1cc62d39a381/plugin?embedAsThumbnail=false&displayTitle=false&startTime=0&autoPlay=false&hideControls=true&showCaptions=false&width=1280&height=720&displaySharing=false&displayNotes=false&displayAttachments=false&displayLinks=false&displayEmbedCode=false&displayDownloadIcon=false&displayMetaData=false&displayCredits=false&displayCaptionSearch=false&audioPreviewImage=false&displayViewersReport=false&displayAxdxs=false&displayComments=false "BeadEvaluator")


## How To Build The Plugin

### 1. Dependencies

- JDK >= 1.8 (tested with JDK 8-13)
- Internet Connection (other dependencies will be downloaded automatically)
- Optional: IDE with [Gradle](http://www.gradle.org/) support


### 2. Configuration (Optional)

If the plugin shall be installed to a custom destination, specify the path in `build.properties`, e.g.,
    
    # vrl property folder location (plugin destination)
    vrldir=/path/to/.vrl/0.4.4/myvrl
    
Otherwise, the plugin will be installed to the default location (depends on VRL version that is specified in the gradle dependencies).

### 3. Build & Install

#### IDE

To build the project from an IDE do the following:

- open the  [Gradle](http://www.gradle.org/) project
- call the `installVRLPlugin` Gradle task to build and install the plugin
- restart VRL-Studio

#### Command Line

Building the project from the command line is also possible.

Navigate to the project folder and call the `installVRLPlugin` [Gradle](http://www.gradle.org/)
task to build and install the plugin.

##### Bash (Linux/OS X/Cygwin/other Unix-like OS)

    cd Path/To/Biofilm
    bash ./gradlew installVRLPlugin
    
##### Windows (CMD)

    cd Path\To\Biofilm
    gradlew installVRLPlugin

Finally, restart VRL-Studio
