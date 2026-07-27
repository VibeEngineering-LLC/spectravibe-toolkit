/*
	spheader.h

	Формат файла спектра создаваемого демонстрационной программой
*/

#ifndef __SPHEADER_H
#define __SPHEADER_H

#pragma pack(1)

typedef struct
{
	BYTE device[ 16 ];		/* имя устойства */
	BYTE unit;				/* номер АЦП */
	BYTE section;	  		/* номер секции */

	BYTE sampleTime[ 6 ];  	/* время взятия образца "ДМГЧМС" */
	BYTE weight[ 16 ];    	/* масса, кг */
	BYTE volume[ 16 ];		/* объем, л */

	BYTE startTime[ 6 ];    	/* время начала набора "ДМГЧМС" */
	DWORD liveTime;			/* живое время */
	DWORD realTime;			/* реальное время */

	WORD buffer;			/* размер буфера АЦП */
	WORD gain;	      		/* разрешение АЦП */
	WORD offset;			/* смещение АЦП */
	WORD lowerLevel;		/* ДНУ */
	WORD upperLevel;		/* ДВУ */

	WORD first;	      		/* первый канал */
	WORD last;	      		/* последний канал */

	WORD chWidth;			/* ширина канала */
	WORD adcCont[ 10 ];		/* продолжение параметров АЦП */

	BYTE prepar;			/* тип подготовки: 0-нет, 1-озол */
	BYTE preparTime[ 6 ];   	/* время начала отбора "ДМГЧМС" */

	BYTE reserve[ 15 ];

	BYTE energy[ 4 ][ 16 ];		/* коэффициенты калибровки по энергии */
	BYTE fwhm[ 4 ][ 16 ];		/* коэффициенты калибровки по полуширине */

	BYTE comment[ 4 ][ 64 ];	/* описание условий измерения */
} tSpectrHeader, *pSpectrHeader;

typedef struct
{
      tSpectrHeader header;
      DWORD *data;
}tSpectr;

#pragma pack()

#endif
