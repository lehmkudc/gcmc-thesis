This collection of runs is to characterize how repeatable these runs are. By "Repeatable" I mean how much equilibration runs tend to relate to eachother. 


# Experiment 1:
* 20 single-runs at P=100 yco=[0, 0.2, 0.4, 0.6, 0.8, 1.0]
* 20 single-runs at P=150 yco=[0, 0.2, 0.4, 0.6, 0.8, 1.0]
* 20 single-runs at P=200 yco=[0, 0.2, 0.4, 0.6, 0.8, 1.0]

Making sure to give each run plenty of room to equilibrate. I will determine the appropriate number of iterations by running each attempt one at a time. 


# Setup:
For each of the 18 p/y
  * have a script that is running on repeat 20 times
  * Determine for each one how repeatable their results are. Do runs tend to equilibrate to different local minima
  * Determine if there is a systematic offshoot of resulting pressures and yco at different presssures and compositions
  * Ideally this would confirm (or give an idea as to how off we were) our Fugacity calculations.
  * If there is a run with not enough equilibration time, I'll be able to notice and set a higher run count

  