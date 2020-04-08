FROM node:12.16.1

WORKDIR /app

COPY . .

RUN yarn install --production=false

EXPOSE 5000

CMD ["yarn", "dev"]
