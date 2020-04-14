# Epidemic Calculator 

## Overview
This is a modified verision of Gabrielle Goh's exellent epidemic calculator found [here](http://gabgoh.github.io/COVID/) ([source code](https://github.com/gabgoh/epcalc)). The first major change is that you can now extend the number of days that the model computes out to >18 months. The second major change is that you can now specify how long the social distancing measures are implemented, followed by how long the social distancing measures are relaxed. This is so that the Covid-19 "aftershocks" can be studied. The default is cycles of social distancing intervention for 90 days (grey shaded regions), followed by 30 days of returning to life as normal. You can also change the R0 for the time periods where social distancing is relaxed. The initial parameters are based roughly on those estimated by the [CDC for Wuhan](https://wwwnc.cdc.gov/eid/article/26/7/20-0282_article), however, one should be careful when trying to draw too strong a conclusion from the model as many parameters are not well known and can fluctuate wildly depending on the time and place of the outbreak.

# Installation  
The easiest install method is to use Docker and docker-compose. Clone this repo onto your computer, then run:

```docker build .```  
```docker-compose up```  

This will start a local server on your machine on port 5000. Point a browser to [localhost:5000](http://localhost:5000) to access it. You can find the code in the ```src``` directory. After editing, the code will automatically recompile, and force refreshing the browser will reflect the changes.

To stop the program, just run:  

```docker-compose down```