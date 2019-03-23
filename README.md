# Hybrid-Storage_Project
Project on the optimization of an electric vehicle's energy storage system using machine learning, to help EV drivers.

## Background
Range anxiety is a major issue that discourages drivers from considering electric cars as an alternative mode of transportation to the standard combustion engine vehicle. This project aims to provide such drivers with an upper bound prediction of travel time if they were to use an electric vehicle based on their own driving habits. The goal is to help such drivers estimate whether it would be reasonable for them to considering to this more environmentally-friendly alternative without impacting their current driving habits.

## Project
In this project, a battery storage system in electric cars is being simulated under a very basic model of energy transfers. (This is a very common system consisting of two storage devices, and can be found in many cars today.) AI is being used in the form of reinforcement learning, to determine the most optimal electronic control of the energy transfers for the given driving route. By learning from the user's historical driving patterns, the software agent is able to both minimize energy losses and increase the battery lifetime, to allow the driver to travel further in the future. As far as we know, such techniques are currently only at the research stage today, but this simulator will approximate the time with some accuracy provided that the driver can input his/her driving speed data (which could be easily measured from an odometer during the trip).

A secondary problem that is also being solved the optimal sizing of the energy storages to maximize range while minimizing cost.

## Theory
The RL framework being used to solve the problem is a planning technique called of approximate dynamic programming (ADP). The ADP has been reformulated as a linear program to allow for fast convergence to an optimal solution. Feature vectors have been developed based on visualizing the optimal value function for a training dataset.

The environment is a very basic physical model of two energy storages in the electric vehicle that treats them as wells, extended from a mathematical queuing model. (Specifically, this a solving a 2-storage inventory control problem.) Although we are ignoring many losses and making some assumptions, the time for fully depleting the stored energy under optimal conditions gives an estimate of the driving range under idealized conditions. In current visualizations of the trip, the state of charge (energy) over time is what is shown.
At present, the driving data is being randomly generated following a beta distribution in the environment and simulator. This will be replaced by a categorical distribution once the training dataset is ready.

## Development
The training environment and simulator have been developed in MATLAB, and are currently running on a 4-core home PC with no GPU support currently needed. The program will be able to use historical data to predict the best way to internally control the energy transfers for driving. It is currently running in real-time at a control period of 0.1s and typically 10,000+ states during the optimization process. We are in the process of porting the code to Python and then parallelizing it for speed.

No common standards or toolkits have used in developing this agent, as this is primarily an academic project with not much other work done using the ADP optimizer. However, we will be using PyTorch to develop an MLP to approximate the value function, as an improvement on the current linear approximation.