#!make
IMAGE_NAME=skasa/hi-inator

ifndef config
    config=parameters.json
endif

.PHONY: all build run force-build

all: build run


build:
	docker build -t $(IMAGE_NAME) .

force-build:
	docker build --pull -t $(IMAGE_NAME) --no-cache=true .

run:
	docker run -v `pwd`/input:/input:ro -v `pwd`/output:/output:rw -e CONFIG=/input/$(config) $(IMAGE_NAME) 

shell:
	docker run -v `pwd`/input:/input:ro -v `pwd`/output:/output:rw -e CONFIG=/input/$(config) -it $(IMAGE_NAME)  bash
