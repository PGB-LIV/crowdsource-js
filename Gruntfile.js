module.exports = function(grunt) {

  // Project configuration.
  grunt.initConfig({
    pkg: grunt.file.readJSON('package.json'),
    clean: ['dist/*.js', 'test/testem.tap'],
	mkdir: {
		all: {
			options: {
				mode: 0700,
				create: ['dist/cov', 'dist/metrics']
			}
		}
	},
    jshint: {
      all: ['src/*.js'],
      options: grunt.file.readJSON('build/jshint.json')
    },
    concat: {
      build: {
        files: {
          'dist/<%= pkg.name %>.js': [
            'src/crowdsource.js',
            'src/cs_worker.js',
            'src/thirdphase_worker.js'
          ]
        }
      }
    },
	closurecompiler: {
		all: {
		  options: {
		    "closure_compilation_level": 'ADVANCED'
		  },
		  files: {
            'dest/closure_compile_adv.min.js': ['src/*.js']
          }
		}		
    },
    uglify: {
      options: {
        banner: '/*! <%= pkg.name %> <%= grunt.template.today("yyyy-mm-dd") %> */\n'
      },
      build: {
        src: 'dist/<%= pkg.name %>.js',
        dest: 'dist/<%= pkg.name %>.min.js'
      }
    },
    testem: {
      options: {
        launch_in_ci: ['PhantomJS']
      },
      'test/testem.tap': ['test/*.html']
    },
    "qunit-cov": {
      test: {
        minimum: 0.9,
        srcDir: 'src',
        depDirs: ['test'],
        outDir: 'dist/cov',
        testFiles: ['test/*.html']
      }
    },
    plato: {
      options: {
        title: 'CrowdSourcingClient-JS',
        jshint: grunt.file.readJSON('build/jshint.json')
      },
      metrics: {
        files: {
          'dist/metrics': [ 'src/*.js' ]
        }
      }
    }
  });

  grunt.loadNpmTasks('grunt-contrib-clean');
  grunt.loadNpmTasks('grunt-mkdir');
  grunt.loadNpmTasks('grunt-contrib-jshint');
  grunt.loadNpmTasks('grunt-contrib-concat');
  grunt.loadNpmTasks('grunt-contrib-uglify');
  grunt.loadNpmTasks('grunt-testem');
  grunt.loadNpmTasks('grunt-qunit-cov');
  grunt.loadNpmTasks('grunt-plato');
  grunt.loadNpmTasks('grunt-google-closure-tools-compiler');

  // Default task(s).
  grunt.registerTask('default', ['jshint', 'testem', 'clean', 'mkdir', 'closurecompiler', 'qunit-cov']);
  grunt.registerTask('jenkins', ['jshint', 'testem', 'clean', 'mkdir', 'closurecompiler', 'qunit-cov', 'plato', 'concat', 'uglify']);

};