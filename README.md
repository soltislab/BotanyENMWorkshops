# BotanyENMWorkshop2026
**Stay tuned for the Botany 2026 update!**

This repository contains material from the Botany 2026 workshop titled 'Using Digitized Collections-Based Data in Research: Applications for Ecology, Phylogenetics, and Biogeography' which was held in August 2026.

HTML version of our R markdown can also be viewed here: https://soltislab.github.io/BotanyENMWorkshops/

## Organizers 
Current organizers: Pam Soltis, Doug Soltis, Shelly Gaynor, Makenzie Mabry, Elizabeth White, JT Miller, Sarah Ellen Strickland, and Cameron McMullen.

## Workshop participants:    
- During the workshop, we will work on HiperGator. No need to clone the repository locally.     
- To get started, navigate to `Botany2026/` and open the .Rproj object. Start with 01_Setup.R. 
- If you join the workshop late or if you are using this material outside of the workshop, please clone this repository locally. Then navigate to `Botany2026/` and open the .Rproj object. Start with 01_Setup.R, here you will also download all the data needed to run these tutorials.

### Data information:    
The data for this workshop has been moved to Zenodo. Please download the data folder from Zenodo: [doi.org/10.5281/zenodo.16755492](https://zenodo.org/records/16755492) 

## rJava issues

### M2 debugging - Shelly
It seems that dismo only likes a certain version of openjdk. 

I ended up having to do the following on my M2 mac which has homebrew installed and is using bash (not zsh):

```
brew --install openjdk@17
export JAVA_HOME="$(brew --prefix openjdk@17)/libexec/openjdk.jdk/Contents/Home"
export PATH="$JAVA_HOME/bin:$PATH"
source ~/.bash_profile
JDK_PATH="$(brew --prefix openjdk@17)/libexec/openjdk.jdk/Contents/Home"
sudo R CMD javareconf -e JAVA_HOME="$JDK_PATH"
```

Then uninstall dismo and rJava (`remove.packages("dismo", "rJava")`), restart R, and reinstall dismo and rJava (`install.packages("rJava", "dismo")`) . 

To verify the java version in R was what I wanted, I did `system("java -version")` in R.
