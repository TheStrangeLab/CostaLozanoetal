Note: A newer, more advanced version is in the works and is
available in the "eeg_ana" CVS project.  I highly recommend checking
that version out, as I am rapidly debugging it and it should become
stable soon.

NWM 
3/9/08

This README gives information on a series of scripts designed to
streamline eeg analysis, in conjunction with other scripts in the
eeg_toolbox.  They are a work in progress.

1. First off, all the standard post-processing should be done.  That
   is, the eeg data should be rereferenced, an events struct should be
   created from the behavioral data, and the data should be aligned.
   If this is a scalp study, bad channel and artifact detection should
   be completed and the necessary files saved; if this is an
   intracranial study, good_leads.txt should be created.

2. Next, you need to create a struct that contains the identifiers for
   all the subjects in your experiment, as well as the locations where
   their events.mat files can be found.  The struct can either be
   saved in a .mat file, or you can write an m-file that creates
   the struct.
   
   For example, the m-file could look like this:
     subj(1).id = 'UP001';
     subj.session(1).eventsFile = '/data/eeg/UP001/catFR/session_0/events.mat';
     subj.session(2).eventsFile = '/data/eeg/UP001/catFR/session_1/events.mat';
     subj(2).id = 'UP002';
   and so forth.

3. Run either init_scalp or init_iEEG to produce the 'eeg' struct that
   will hold all the necessary information about the experiment.  This
   struct will be saved in whatever 'resDir' you specify.  Of note is
   the 'params' argument for the 'init_' scripts.  Params is a struct
   that holds all the parameters used for your analysis, which is
   saved as eeg.params, so it is easy to go back and double-check how
   results were generated.  The basic fields of params are:

   eventFilter - a string that is passed into filterStruct.  each
   event struct used for analysis is first filtered with this string.
   This will depend on which events you are interested in.

   binSizeMS - allows you to average the data over a certain size of
   time bin.
   relativeMS - contains the beginning and end of the baseline period
   for each event; used for baseline subtraction and z-transforming
   bufferMS - time window before and after an event to use for
   filtering
   filtfreq - range of frequencies to filter
   filttype - type of filter to use
   filtorder - filter order
   resampledRate - rate to resample data at

4. pass 'eeg' into one of the following scripts:

   create_volt_pattern.m
   create_pow_pattern.m

   This will save out a .mat file for each subject containing either
   voltage or power values for all time and/or frequency bins for each
   subject.  The path to each file will be saved in eeg.subj(n).patFile.

5. pass 'eeg' and 'regParams' into make_regressors.m.  This will
   create regressors for each subject, using whatever 'events'
   fieldnames you specify in regParams.  The path to each file will be
   saved in eeg.subj(n).regFile.

6. Now that you have patterns and regressors saved, you are ready to
   analyze!  Try these:

   pow_anova.m (in progress)

   Or write your own!  It's fun!

Misc Functions

loadStruct.m
Used to change all filenames contained within the eeg struct.  Useful
for moving from hippo to a local machine.

structDefaults.m
Used to set default values for fields in the params struct.
   
