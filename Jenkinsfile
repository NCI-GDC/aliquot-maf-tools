#!groovy

PROJECT_NAME = "aliquot-maf-tools"
PYENV_VERSION = "3.6-dev"

BRANCH_NAME = env.BRANCH_NAME
GIT_HASH = ""
VERSION = ""

PROXY = "http://cloud-proxy:3128"

pipeline {
  agent any
  environment {
        TWINE_REPOSITORY_URL = credentials("${BRANCH_NAME == 'main' ? 'twine_repository_url_prod' : 'twine_repository_url'}")
	TWINE_USERNAME = credentials('twine_username')
	TWINE_PASSWORD = credentials('twine_password')
	QUAY_USERNAME = credentials('QUAY_USERNAME')
	QUAY_PASSWORD = credentials('QUAY_PASSWORD')
  }
  options {
    disableConcurrentBuilds()
    skipStagesAfterUnstable()
  }

  stages {
    stage('Docker Build') {
      steps {
        vbash "make build-docker PROXY=${PROXY}"
      }
    }
    stage('Docker Test') {
      steps {
        sh 'make test-docker'
      }
    }
    stage('Docker Publish') {
      steps {
        script {
          DOCKER_IMAGE = sh(script: "make version-docker-tag", returnStdout: true).trim()
	}
        sh "make publish-docker DOCKER_IMAGE=${DOCKER_IMAGE}"
      }
    }
  }
  post {
    always {
      sh 'make clean'
    }
  }
}

def vbash(command) {
  sh """#!/bin/bash
        eval \"\$(pyenv init -)\"
	eval \"\$(pyenv virtualenv-init -)\"

	pyenv virtualenv ${PYENV_VERSION} ${PROJECT_NAME}-venv || true
	pyenv activate ${PROJECT_NAME}-venv

	${command}
  """
}
