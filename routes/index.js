var express = require('express');
var router = express.Router();

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('test', {
    title: 'Red Nova',
    hasResults: false
  });
});

/* POST compute page. */
router.post('/compute', function(req, res, next) {

  try {

    if (!req.body.data) {
      res.redirect('/');
    }

    var options = {
      mode: 'text',
      pythonPath: 'python',
      pythonOptions: ['-u'],
      scriptPath: '',
      args: ["KIC 9832227" ,req.body.data]
    };

    console.log(req.body.data);

    PythonShell.run('OrbitalPhasePackage.py', options, function (err, OrbData) {
      // if (err) throw err;
      if (err) {
        res.json({
          errorType:"Python Error",
          error: err
        });
      };
      // results is an array consisting of messages collected during execution
      // res.json(OrbData);

      // res.json({
      //   title: 'Red Nova',
      //   hasResults: true,
      //   results:OrbData
      // });

      res.render('test', {
        title: 'Red Nova',
        hasResults: true,
        results:OrbData
      });

    });
  }
  catch (err) {
    res.json({error:"error"});
  }
  //
  // res.json({
  //   test:req.body.name
  // });
});

module.exports = router;
