FROM penthesilea/imager

RUN pip install psutil simms

# Remove native recipe
RUN rm /code -rf

ADD src /code

RUN mkdir -p /input /output

ENV INPUT /input
ENV OUTPUT /output

WORKDIR /code
CMD  sh run.sh
