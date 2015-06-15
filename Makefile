#!make
IMAGE_NAME=skasa/hi-inator

ifndef config
    config=parameters.json
endif

.PHONY: all build run force-build

all: download build run

download:
	./download.sh

build:
	./download.sh && docker build -t $(IMAGE_NAME) .

force-build:
	./download.sh && docker build --pull -t $(IMAGE_NAME) --no-cache=true .

run:
	docker run -v `pwd`/input:/input:ro -v `pwd`/output:/output:rw -e config=$(config) $(IMAGE_NAME) 

shell:
	docker run -ti -v `pwd`/input:/input:ro -v `pwd`/output:/output:rw $(IMAGE_NAME) bash .
