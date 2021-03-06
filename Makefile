REPO = aliquot-maf-tools

MODULE = aliquotmaf

COMMIT_HASH:=$(shell git rev-parse HEAD 2>/dev/null)

DOCKER_REPO := quay.io/ncigdc
DOCKER_IMAGE_COMMIT := ${DOCKER_REPO}/${REPO}:${COMMIT_HASH}
DOCKER_IMAGE_LATEST := ${DOCKER_REPO}/${REPO}:latest

TWINE_REPOSITORY_URL?=
PIP_EXTRA_INDEX_URL?=

.PHONY: version version-* print-*
version:
	@echo --- VERSION: ${PYPI_VERSION} ---

version-docker:
	@python setup.py -q print_version --docker

version-docker-tag:
	@docker run --rm --entrypoint="make" ${DOCKER_IMAGE_LATEST} "version-docker"

.PHONY: docker-login
docker-login:
	docker login -u="${QUAY_USERNAME}" -p="${QUAY_PASSWORD}" quay.io


.PHONY: build build-* clean clean-* init init-* lint requirements run version
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

clean-docker:
	@echo "Removing latest image: ${DOCKER_IMAGE_LATEST}"
	docker rmi -f ${DOCKER_IMAGE_LATEST}

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
	docker build . \
		--file ./Dockerfile.multistage \
		--build-arg http_proxy=${PROXY} \
		--build-arg https_proxy=${PROXY} \
		-t "${DOCKER_IMAGE_COMMIT}" \
		-t "${DOCKER_IMAGE_LATEST}"

build-pypi:
	@echo
	@echo Building wheel - ${PYPI_VERSION}
	python3 setup.py bdist_wheel -b ${MODULE}.egg-info

.PHONY: test test-* tox
test: lint test-unit

tox:
	@tox

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

.PHONY: publish-*
publish-docker:
	docker tag ${DOCKER_IMAGE_COMMIT} ${DOCKER_REPO}/${REPO}:${DOCKER_TAG}
	docker push ${DOCKER_IMAGE_COMMIT}
	docker push ${DOCKER_REPO}/${REPO}:${DOCKER_TAG}

publish-pypi:
	@echo
	@echo Publishing wheel
	python3 -m twine upload $(shell ls -1 dist/*.whl | head -1)
