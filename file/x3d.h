/**
 * simplified x3d file handling
 * @author tweber
 * @date nov-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __X3D_FILES_H__
#define __X3D_FILES_H__

#include <ostream>
#include <vector>

#include "../math/linalg.h"
#include "../math/quat.h"

namespace tl{


class X3dElem;
class X3dElem
{
	protected:
		std::vector<X3dElem*> m_vecChildren;

	public:
		X3dElem() = default;
		virtual ~X3dElem();

		virtual void Write(std::ostream& ostr) const;
		virtual void AddChild(X3dElem* pElem);
};

class X3dScene : public X3dElem
{
	public:
		X3dScene() = default;
		virtual ~X3dScene() = default;

		virtual void Write(std::ostream& ostr) const;
};

class X3dTrafo : public X3dElem
{
	public:
		typedef double t_real;
		typedef ublas::vector<t_real> t_vec;
		typedef ublas::matrix<t_real> t_mat;
		typedef math::quaternion<t_real> t_quat;

	protected:
		t_vec m_vecTrans;
		t_vec m_vecScale;
		t_quat m_quatRot;
		bool m_bHasRot = 0;

	public:
		X3dTrafo() = default;
		virtual ~X3dTrafo() = default;

		virtual void Write(std::ostream& ostr) const override;

		void SetTrans(const t_vec& vec) { m_vecTrans = vec; }
		void SetScale(const t_vec& vec) { m_vecScale = vec; }
		void SetRot(const t_quat& quat) { m_quatRot = quat; m_bHasRot = 1; }
};

class X3dSphere : public X3dElem
{
	public:
		typedef double t_real;
		typedef ublas::vector<t_real> t_vec;

	protected:
		t_real m_dRadius = 1.;
		t_vec m_vecColor;

	public:
		X3dSphere() = default;
		X3dSphere(t_real dRad) : m_dRadius(dRad) {}
		virtual ~X3dSphere() = default;

		virtual void Write(std::ostream& ostr) const override;

		void SetRadius(t_real dRad) { m_dRadius = dRad; }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }
};

// -----------------------------------------------------------------------------

class X3d
{
	protected:
		X3dScene m_scene;

	public:
		X3d() = default;
		virtual ~X3d();

		void Write(std::ostream& ostr) const;
		bool Save(const char* pcFile) const;

		X3dScene& GetScene() { return m_scene; }
		const X3dScene& GetScene() const { return m_scene; }
};

}
#endif
