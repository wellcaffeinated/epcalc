jquery('#interactive-byline .last-byline').after('<div>Model created with Gabriel Goh, Steven De Keninck, Ashleigh Tuite and David N. Fisman</div>');

// Resize...
jquery('.model-wrapper.single .model-wrap').css({
  width:jquery('.model-wrapper.full .model-wrap').outerWidth()
});
var keys = {
  'main': ['Dead','Hospital','Infected']//,'Exposed']
};
var overallN = 150;

function makeViz(){

  var margin = {
    top: 100,
    left: 15,
    bottom: 40,
    right: 75
  };

  var height =  window.innerWidth < 600 ? window.innerHeight*.6 - margin.top - margin.bottom : window.innerWidth < 900 ? window.innerHeight * .8 - margin.top - margin.bottom :   window.innerHeight*.8 - margin.top - margin.bottom;
  var width = jquery('.model-wrapper.full .model-viz').outerWidth() - margin.left-margin.right;

  // Main SVG
  var svgW = select('.model-wrapper.full .model-viz').append('svg')
    .attr('width', width+margin.left+margin.right)
    .attr('height', height+margin.top+margin.bottom);
  var svg = svgW.append('g')
    .attr('transform', `translate(${margin.left}, ${margin.top})`);

  // Lines and scales
  var startDate = new Date('2020/01/01');
  var endDate = new Date('2020/08/06');
  var xDate = time().domain([startDate, endDate]).range([0, width]);
  var x = linear$1().domain([0, 220]).range([0, width]);
  // var y = d3.scaleLinear().range([0, height])
  var yInvert = linear$1().range([height,0]);
  // var opacity = d3.scaleLinear().domain([height*.8,height*.9]).range([1,0])
  // var textPath = d3.line()
  //   .x(d => xDate(new Date(d.Time)))
  //   .y(d => yInvert(d.y))
  // var line = d3.line()
  //   .x(d => xDate(new Date(d.date)))
  //   .y(d => yInvert(d.cumulative_cases*caseMultiplier))

  var mos = 'Jan.,Feb.,March,April,May,June,July,Aug.,Sept.,Oct.,Nov.,Dec.'.split(',');
  var xFormat = (d,i) => {
    return mos[moment(d).format('M')-1]
  };
  var yFormat = (d,i) => {
    return d > 999999999 ? (d/1000000000).toFixed(1).toLocaleString() + 'B' : d > 999999 ? (d/1000000).toFixed(1).toLocaleString() + 'M' : d > 999 ? (d/1000).toFixed(1).toLocaleString() + 'K' : d
  };
  var totalYFormat = (d,i) => {
    var n = d > 999999999 ? ((Math.round((d/100))*100)/1000000000).toFixed(1) + ' billion' : d > 999999 ? ((Math.round((d/100))*100)/1000000).toFixed(1) + ' million' : d > 999 ? (((Math.round((d/100))*100)/1000)*1000).toLocaleString() : Math.round(d);
    return n.replace('.0','')
  };

  // Axes
  var xAxis = svg.append('g')
    .attr('transform', `translate(0, ${height})`)
    .call(axisBottom(x).tickFormat(xFormat).ticks(4))
    .attr('class', 'x-axis');
  var yAxisCreator = axisRight(yInvert).tickFormat(yFormat).tickSize(width).ticks(5);
  var yAxis = svg.append('g')
    .call(yAxisCreator)
    .attr('class', 'y-axis');



/* - DATA AND THE MODEL - */
  var data = [];
  var stack$1 = [];
  var totals = [];

  var usPop = 327000000*.7;
  // Defaults
  var R0                = 2.5;
  var CFR               = 0.01;
  var InterventionAmt   = .4;
  var InterventionTime  = 73;
  var P_SEVERE          = 0.2;
  var N                 = usPop;///6999999.999999998 //Math.exp(logN)
  var duration          = 14;
  var seasonal_effect   = .46;
  var old_seasonal_effect = .46;
  var location  = 'United States'; // world
  function processData(){
    if(N > usPop) seasonal_effect = 0;

    function f() {

      var Time_to_death     = 32;
      // var N                 = 327000000
      var I0                = {
        'United States': 3.4,
        'world': 107
      };
        // - (10,048 - 3248) = 6800 deaths on March 19th // USA - 3.4  - 219 deaths on March 19th
      // var R0                = 2.2
      var D_incbation       = 5.2;
      var D_infectious      = 2.9;
      var D_recovery_mild   = (14 - 2.9);
      var D_recovery_severe = (31.5 - 2.9);
      var D_hospital_lag    = 5;
      var D_death           = Time_to_death - D_infectious;
      var dt                = 2;
      // var P_SEVERE          = 0.2
      // var duration          = 7*12*1e10
      // var seasonal_effect          = 1

      var interpolation_steps = 4;
      var steps = overallN*interpolation_steps;
      var dt = dt/interpolation_steps;
      var sample_step = interpolation_steps;
      Date.prototype.addDays = function(days) {
        var date = new Date(this.valueOf());
        date.setDate(date.getDate() + days);
        return date;
      };
      Date.prototype.isLeapYear = function() {
          var year = this.getFullYear();
          if((year & 3) != 0) return false;
          return ((year % 100) != 0 || (year % 400) == 0);
      };

      // Get Day of Year
      Date.prototype.getDOY = function() {
          var dayCount = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
          var mn = this.getMonth();
          var dn = this.getDate();
          var dayOfYear = dayCount[mn] + dn;
          if(mn > 1 && this.isLeapYear()) dayOfYear++;
          return dayOfYear;
      };

      function dateFromDay(t) {
        var date = new Date("01/01/2020"); // Initial date
        return date.addDays(t)
      }

      // f is a func of time t and state y
      // y is the initial state, t is the time, h is the timestep
      // updated y is returned.
      var integrate=(m,f,y,t,h)=>{
        for (var k=[],ki=0; ki<m.length; ki++) {
          var _y=y.slice(), dt=ki?((m[ki-1][0])*h):0;
          for (var l=0; l<_y.length; l++) for (var j=1; j<=ki; j++) _y[l]=_y[l]+h*(m[ki-1][j])*(k[ki-1][l]);
          k[ki]=f(t+dt,_y,dt);
        }
        for (var r=y.slice(),l=0; l<_y.length; l++) for (var j=0; j<k.length; j++) r[l]=r[l]+h*(k[j][l])*(m[ki-1][j]);
        return r;
      };

      var Integrators = {
        Euler    : [[1]],
        Midpoint : [[.5,.5],[0, 1]],
        Heun     : [[1, 1],[.5,.5]],
        Ralston  : [[2/3,2/3],[.25,.75]],
        K3       : [[.5,.5],[1,-1,2],[1/6,2/3,1/6]],
        SSP33    : [[1,1],[.5,.25,.25],[1/6,1/6,2/3]],
        SSP43    : [[.5,.5],[1,.5,.5],[.5,1/6,1/6,1/6],[1/6,1/6,1/6,1/2]],
        RK4      : [[.5,.5],[.5,0,.5],[1,0,0,1],[1/6,1/3,1/3,1/6]],
        RK38     : [[1/3,1/3],[2/3,-1/3,1],[1,1,-1,1],[1/8,3/8,3/8,1/8]]
      };

      var method = Integrators["RK4"];

      var offset = - dateFromDay(0).getDOY() + (new Date("02/01/2020")).getDOY(); // Coldest day of year

      function f(t, x){

        // SEIR ODE
        // SEIR ODE
        if (t > InterventionTime && t < InterventionTime + duration){
          var beta = (InterventionAmt)*R0/(D_infectious);
        } else if (t > InterventionTime + duration) {
          var beta = 1*R0/(D_infectious);
        } else {
          var beta = R0/(D_infectious);
        }

        var forcing = (t) => (1 + seasonal_effect*Math.cos(2*3.14159265*(Math.floor(t) - offset)/365))
        beta = beta*forcing(t)/forcing(0) // Forcing, with R0 correction

        var a     = 1/D_incbation;
        var gamma = 1/D_infectious;

        var S        = x[0]; // Susectable
        var E        = x[1]; // Exposed
        var I        = x[2]; // Infectious
        var Mild     = x[3]; // Recovering (Mild)
        var Severe   = x[4]; // Recovering (Severe at home)
        var Severe_H = x[5]; // Recovering (Severe in hospital)
        var Fatal    = x[6]; // Recovering (Fatal)
        var R_Mild   = x[7]; // Recovered
        var R_Severe = x[8]; // Recovered
        var R_Fatal  = x[9]; // Dead

        var p_severe = P_SEVERE;
        var p_fatal  = CFR;
        var p_mild   = 1 - P_SEVERE - CFR;

        var dS        = -beta*I*S;
        var dE        =  beta*I*S - a*E;
        var dI        =  a*E - gamma*I;
        var dMild     =  p_mild*gamma*I   - (1/D_recovery_mild)*Mild;
        var dSevere   =  p_severe*gamma*I - (1/D_hospital_lag)*Severe;
        var dSevere_H =  (1/D_hospital_lag)*Severe - (1/D_recovery_severe)*Severe_H;
        var dFatal    =  p_fatal*gamma*I  - (1/D_death)*Fatal;
        var dR_Mild   =  (1/D_recovery_mild)*Mild;
        var dR_Severe =  (1/D_recovery_severe)*Severe_H;
        var dR_Fatal  =  (1/D_death)*Fatal;

        //      0   1   2   3      4        5          6       7        8          9
        return [dS, dE, dI, dMild, dSevere, dSevere_H, dFatal, dR_Mild, dR_Severe, dR_Fatal]
      }

      var v = [1, 0, I0[location]/(N-I0[location]), 0, 0, 0, 0, 0, 0, 0];
      var t = 0;

      var P     = [];
      while (steps--) {
        if ((steps+1) % (sample_step) == 0) {
          // Stuart moves
                P.push({"Time":dateFromDay(t),
                        "Dead":N*v[9],
                        "Susceptible": N*(v[0]),
                        "Hospital": N*(v[5]+v[6]),
                        "Recovered": N*(v[7] + v[8]),
                        "Infected": N*(v[2]+v[3]+v[4]),
                        "Exposed": 0,///N*v[1],
                        "Mild": N*v[3],
                        "Severe": N*v[4],
                        "Severe_H": N*v[5],
                        "Fatal": N*v[6],
                        "Sum": N*v.reduce((a, b) => a + b, 0)});
        }
        v =integrate(method,f,v,t,dt);
        t+=dt;
      }
      return P
    }

    // * - stuart here - */

    data = f();
    data.map(d => {
      d.tot = d.Dead + d.Hospital + d.Infected;// + d.Exposed
      d.tothospital = d.Hospital;
      d.totseverity = d.Dead + d.Hospital + d.Infected;// + d.Exposed
      d.totdead = d.Dead;
      return d
    });

    var lastCol = data[data.length-1];
    totals = {
      'infections': lastCol.Sum-lastCol.Susceptible,
      'dead': max(data, d => d.Dead),
      'recoveries': max(data, d => d.Recovered),
      'infections-peak': max(data, d => (d.Infected+d.Hospital)), //d.Exposed+
      'infections-peak-date': data.filter(d => {
        var max$1 = max(data, d => (d.Infected+d.Hospital));//+d.Exposed
        return (d.Infected + d.Hospital) == max$1 //d.Exposed +
      })[0],
      'location': location
    };
    // console.log({totals})
    // console.log({data})

    stack$1 = stack().keys(keys['main'])(data)
      .map((d,i) => {
        return d.map(d=> {
          d.key = keys.main[i];
          return d
        })
      });
  }