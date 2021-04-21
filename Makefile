# UPDATE ME
REPO = aliquot-maf-tools

# UPDATE ME
MODULE = aliquotmaf

GIT_SHORT_HASH:=$(shell git rev-parse --short HEAD)

PYPI_VERSION:=$(shell python3 setup.py -q print_version --pypi)
DOCKER_VERSION:=$(shell python3 setup.py -q print_version --docker)
COMMIT_HASH:=$(shell python3 setup.py -q print_version --hash)

DOCKER_REPO := quay.io/ncigdc
DOCKER_IMAGE := ${DOCKER_REPO}/${REPO}:${DOCKER_VERSION}
DOCKER_IMAGE_COMMIT := ${DOCKER_REPO}/${REPO}:${COMMIT_HASH}
DOCKER_IMAGE_LATEST := ${DOCKER_REPO}/${REPO}:latest

TWINE_REPOSITORY_URL?=""
PIP_EXTRA_INDEX_URL?=

.PHONY: version version-* print-*
version:
	@echo --- VERSION: ${PYPI_VERSION} ---

version-docker:
	@echo ${DOCKER_IMAGE}
	@echo ${DOCKER_IMAGE_COMMIT}

.PHONY: docker-login
docker-login:
	docker login -u="${QUAY_USERNAME}" -p="${QUAY_PASSWORD}" quay.io


.PHONY: build build-* clean init init-* lint requirements run version
init: init-pip init-hooks

# Include next line if publishing to Jenkins
# --extra-index-url ${TWINE_REPOSITORY_URL}
init-pip:
	@echo
	@echo -- Installing pip packages --
	pip3 install \
		--no-cache-dir \
		-r dev-requirements.txt \
		-r requirements.txt
	python3 setup.py develop

init-hooks:
	@echo
	@echo -- Installing Precommit Hooks --
	pre-commit install

init-venv:
	@echo
	PIP_REQUIRE_VIRTUALENV=true pip3 install --upgrade pip-tools

clean:
	rm -rf ./build/
	rm -rf ./dist/
	rm -rf ./*.egg-info/

lint:
	@echo
	@echo -- Lint --
	python3 -m flake8 ${MODULE}/

run:
	bin/run

requirements: init-venv requirements-prod requirements-dev

requirements-dev:
	python3 setup.py -q capture_requirements --dev
	pip-compile -o dev-requirements.txt dev-requirements.in

requirements-prod:
	pip-compile -o requirements.txt

.PHONY: build build-*

build: build-docker

build-docker:
	@echo
	@echo -- Building docker --
	python3 setup.py build
	mkdir -p dist
	cp -r build/lib/* dist/
	cp -r bin/ dist/
	cp -f Makefile requirements.txt dev-requirements.txt README.md setup.py dist/
	docker build . \
		--file ./Dockerfile \
		--build-arg http_proxy=${PROXY} \
		--build-arg https_proxy=${PROXY} \
		-t "${DOCKER_IMAGE_COMMIT}" \
		-t "${DOCKER_IMAGE}" \
		-t "${DOCKER_IMAGE_LATEST}"

build-pypi:
	@echo
	@echo Building wheel - ${PYPI_VERSION}
	python3 setup.py bdist_wheel -b ${MODULE}.egg-info

.PHONY: test test-*
test: lint test-unit

test-unit:
	@echo
	@echo -- Unit Test --
	python3 -m pytest --cov-report term-missing \
		--junitxml=build/unit-test.xml \
		--cov=${MODULE} \
		tests/

test-docker:
	@echo
	@echo -- Running Docker Test --
	docker run --rm ${DOCKER_IMAGE_LATEST} test

.PHONY: publish publish-*
publish:
	docker push ${DOCKER_IMAGE_COMMIT}
	docker push ${DOCKER_IMAGE}

publish-pypi: dist/*.whl
	@echo
	@echo Publishing wheel
	python3 -m twine upload $(shell ls -1 dist/*.whl | head -1)
