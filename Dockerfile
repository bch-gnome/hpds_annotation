FROM python:3

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY transform_csq.v3.py /

CMD [ "python", "/transform_csq.v3.py" ]
