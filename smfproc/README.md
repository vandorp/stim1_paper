SMFPROC
=======

View, analyze, and export single-molecule FRET movies, traces and model fits.

- Abstract donor and acceptor timeseries traces from a movie
- Preselection based on a simple noise threshold criterium, to get rid of donor-only 
traces (use direct acceptor excitation at the end of the trace to identify zero-fret
molecules).
- Find optimal HMM fit for each trace using SMART routines.
- Use various criteria based on discrete HMM states to identify 'good' FRET traces


## Version
For use with Python 3. Cross-platform in principle. Developed in Ubuntu 16 using 
Python 3.5. Somewhat tested on Windows 10, basically untested on Mac.

Supplement to
"Conformational dynamics of auto-inhibition in the ER calcium sensor STIM1"
Stijn van Dorp, Ruoyi Qiu, Ucheor B. Choi, Minnie M. Wu, Michelle Yen, 
Michael Kirmiz, Axel T. Brunger, Richard S. Lewis
https://doi.org/10.1101/2020.12.17.423361

## Installation

**WINDOWS**

* Install Anaconda. Full Anaconda contains all the necessary Python modules.
* Copy the smfproc program files to a folder named 'smfproc'
* Add the path containing the 'smfproc' folder to the PYTHONPATH environment 
  variable: Computer > Properties > Environment variables
* Open an Anaconda command window and issue 'python smfrun.py'


**MAC/LINUX**

* Copy the smfproc program files to a folder named 'smfproc'
* Add the path containing the 'smfproc' folder to the PYTHONPATH environment 
  variable:
	- Create/open a file .bash_profile (.profile on MAC) in your home directory.
	- Add the following line to the file:
		export PYTHONPATH=$PYTHONPATH:[folder]
		(E.g. when 'smfproc' is on your desktop, [folder] is ~/Desktop)
	- Refresh the environment, or just open a new terminal window
	- Issue 'echo $PYTHONPATH' to check that the folder has been added
* Open a terminal window and issue 'python smfrun.py'


## Program
The movie images from donor and acceptor channels need to be brought into register
before analysis. A registration file can be constructed using the GUI (see below).
This file can be used for subsequent analyses as long as the camera images haven't
shifted. The file needs to be located in any folder upstream of the movie.

Analysis consists of five steps:

1. laser - obtain an average donor fluorescence intensity profile from all movies 
in a folder to estimate the shape of the laser spot.
2. movies - abstracts donor and acceptor traces from each movie.
3. traces - crosstalk correction, calculate FRET ratio, separate signal from noise,
fit HMM to obtain dicrete states.
4. states - select and categorize traces based on various criteria.
5. output - store data and make histograms.


## Usage
Issuing 'python smfrun.py' starts the GUI. See below for usage.

'python smfrun.py --path [folder]'
Starts the GUI and loads the movies in the specified folder.

'python smfrun.py --path [folder] --analyze'
Full analysis of the movies in specified folder, without opening the GUI.

'python smfrun.py --path [folder] --analyze [step]'
Partial analysis of the movies in specified folder, starting from [step].
[step] can be any of the five options above. E.g. if you want to change some trace
selection parameters, but don't want to redo all the trace fitting (the most time-
consuming part) use --analyze states.

'python smfrun.py --help'
Print all command line options.

Descriptions of the analysis settings are provided as comments in the smfrun.py script.


## GUI
**trace selection/deselection**
Double-click on a trace in the listbox (or use the spacebar) to manually select or 
deselect the trace for analysis.
